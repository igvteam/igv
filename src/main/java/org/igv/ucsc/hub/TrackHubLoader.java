package org.igv.ucsc.hub;

import org.igv.feature.genome.Genome;
import org.igv.feature.genome.GenomeDownloadUtils;
import org.igv.feature.genome.GenomeManager;
import org.igv.feature.genome.HostedGenomes;
import org.igv.feature.genome.load.GenomeConfig;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.ui.IGV;
import org.igv.ui.IGVMenuBar;
import org.igv.ui.action.SelectHubTracksAction;
import org.igv.ui.genome.GenomeListItem;
import org.igv.ui.genome.GenomeListManager;
import org.igv.ui.genome.GenomeSelectionDialog;
import org.igv.ui.util.MessageUtils;
import org.igv.ui.util.UIUtilities;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Loads a track hub entered by URL.  A hub references one or more genomes, each either a UCSC hosted genome
 * (including GenArk) referenced by id, or a custom assembly defined by the hub with a "twoBitPath".  A single hub
 * can mix the two.
 * <p>
 * Single genome hub: if the hub references the current genome its track selection dialog is shown, otherwise IGV
 * switches to the hub's genome -- creating it from the hub if it has a twoBitPath, else loading it as if selected
 * from "Load Genome from Host" and then showing the hub's track selection dialog.  A genome switch is expected
 * from the Genomes menu but not from the Tracks menu, so the latter confirms it with the user first.
 * <p>
 * Multiple genome hub: the user is offered the referenced genomes not already in the genome list, and those selected
 * are added, with default annotation selection deferred to their first load.  Then, if the current genome is one of
 * the hub's genomes its track selection dialog is shown, otherwise the user is offered a switch to one of them,
 * followed by the hub's track selection dialog if the chosen genome is a hosted one.
 * <p>
 * The hub is recorded in the HubRegistry (the user managed hubs) for the genomes it references but does not define,
 * so it appears in the Tracks menu for each of them.  A genome the hub defines with a twoBitPath records the hub in
 * its own genome json instead, so an assembly hub is never a managed hub.  Long running; must not be called on the
 * event thread.
 */
public class TrackHubLoader {

    private static Logger log = LogManager.getLogger(TrackHubLoader.class);

    /**
     * @param url                 hub.txt url
     * @param confirmGenomeSwitch if true, and the hub references a single genome other than the current one, ask
     *                            the user to confirm the genome switch before proceeding.  Loading from the Tracks
     *                            menu passes true; a switch is the expectation from the Genomes menu.
     */
    public static void loadHub(String url, boolean confirmGenomeSwitch) throws IOException {

        Hub hub = HubParser.loadHub(url);
        List<String> genomeIds = hub.getGenomeIds();
        if (genomeIds.isEmpty()) {
            throw new IOException("Hub does not reference any genomes: " + url);
        }
        if (genomeIds.size() == 1) {
            loadSingleGenomeHub(hub, genomeIds.get(0), confirmGenomeSwitch);
        } else {
            loadMultiGenomeHub(hub, genomeIds);
        }
    }

    private static void loadSingleGenomeHub(Hub hub, String genomeId, boolean confirmGenomeSwitch) throws IOException {

        Genome currentGenome = GenomeManager.getInstance().getCurrentGenome();

        if (currentGenome != null && hub.referencesGenome(currentGenome.getUCSCId())) {
            // The genome's own hub (the one it was created from) is recorded in its json, not as a managed hub.
            if (!isGenomeHub(currentGenome, hub)) {
                HubRegistry.addUserHub(hub);
            }
            IGVMenuBar.getInstance().updateMenus(currentGenome);
            selectTracks(hub, currentGenome);
            return;
        }

        // The hub references a different genome.  Switching from the Tracks menu is unexpected, so confirm first.
        if (confirmGenomeSwitch) {
            String current = currentGenome == null ? "none" : currentGenome.getDisplayName();
            boolean proceed = MessageUtils.confirm("This hub references the genome \"" + genomeId +
                    "\", not the current genome (" + current + ").  Switching genomes will unload all tracks.  " +
                    "Switch to " + genomeId + "?");
            if (!proceed) {
                addHubWithoutSwitching(hub, genomeId);
                IGVMenuBar.getInstance().updateMenus(currentGenome);
                return;
            }
        }

        // Either way the genome's own default annotation prompt is skipped and the hub's track selection dialog is
        // opened instead, so the user sees one chooser rather than two.
        if (hub.hasTwoBit(genomeId)) {
            // The hub defines the assembly.  Put it in the genome list, as a declined switch does, and load it.
            GenomeListItem item = addAssemblyGenome(hub, genomeId);
            GenomeManager.getInstance().loadGenomeById(item.getId(), false, true);
        } else {
            // A UCSC hosted genome referenced by id.
            GenomeListItem item = findKnownGenome(genomeId);
            if (item == null) {
                MessageUtils.showMessage("The genome referenced by this hub (" + genomeId + ") is not known to IGV.");
                return;
            }
            loadKnownGenome(item);
            HubRegistry.addUserHub(hub);
            // The load rebuilt the menus before the hub was registered
            IGVMenuBar.getInstance().updateMenus(GenomeManager.getInstance().getCurrentGenome());
        }
        selectTracksForCurrentGenome(hub);
    }

    /**
     * Declining the genome switch keeps the current genome, but the hub is still added.  An assembly hub is its
     * genome's own hub, so it is added by putting that genome in the genome list, as the multi-genome path does
     * (annotations are selected on first load); registering it instead would grey it out in the public hub chooser
     * with no genome to load.  Any other hub is registered, and appears in the Tracks menu when its genome is loaded.
     */
    private static void addHubWithoutSwitching(Hub hub, String genomeId) throws IOException {
        if (hub.hasTwoBit(genomeId)) {
            addAssemblyGenome(hub, genomeId);
        } else {
            HubRegistry.addUserHub(hub);
        }
    }

    /**
     * Put an assembly the hub defines in the genome list, if it is not there yet, and return its list item.  The
     * config is copied: the hub's own is cached for the session and must not carry this load's edits.
     */
    private static GenomeListItem addAssemblyGenome(Hub hub, String genomeId) throws IOException {
        GenomeListItem item = findInGenomeList(GenomeListManager.getInstance().getGenomeItemMap(), genomeId);
        if (item == null) {
            GenomeConfig config = hub.getGenomeConfig(genomeId).copy();
            item = addGenome(hub, new HubGenomeSelectionDialog.Candidate(genomeId, config.getName(), config, null));
        }
        return item;
    }

    private static void loadMultiGenomeHub(Hub hub, List<String> genomeIds) throws IOException {

        GenomeListManager genomeListManager = GenomeListManager.getInstance();
        Map<String, GenomeListItem> genomeList = genomeListManager.getGenomeItemMap();

        // Referenced genomes not already in the genome list, and ids IGV does not know.
        List<HubGenomeSelectionDialog.Candidate> candidates = new ArrayList<>();
        List<String> unknownIds = new ArrayList<>();
        for (String genomeId : genomeIds) {
            if (findInGenomeList(genomeList, genomeId) != null) {
                continue;
            }
            if (hub.hasTwoBit(genomeId)) {
                GenomeConfig config = hub.getGenomeConfig(genomeId).copy();
                candidates.add(new HubGenomeSelectionDialog.Candidate(genomeId, config.getName(), config, null));
            } else {
                GenomeListItem item = findHostedGenome(genomeId);
                if (item == null) {
                    unknownIds.add(genomeId);
                } else {
                    candidates.add(new HubGenomeSelectionDialog.Candidate(item.getId(), item.getDisplayableName(), null, item));
                }
            }
        }

        if (candidates.isEmpty() && !unknownIds.isEmpty()) {
            // Nothing to offer; just report the ids IGV does not know and carry on with the genomes it does.
            MessageUtils.showMessage("Genomes referenced by this hub that are not known to IGV and cannot be added: " +
                    String.join(", ", unknownIds));
        } else if (!candidates.isEmpty()) {
            final HubGenomeSelectionDialog[] holder = new HubGenomeSelectionDialog[1];
            UIUtilities.invokeAndWaitOnEventThread(() -> {
                holder[0] = new HubGenomeSelectionDialog(IGV.getInstance().getMainFrame(), hub, candidates, unknownIds);
                holder[0].setVisible(true);
            });
            HubGenomeSelectionDialog dlg = holder[0];
            if (dlg.isCanceled()) {
                return;
            }
            for (HubGenomeSelectionDialog.Candidate candidate : dlg.getSelectedCandidates()) {
                addGenome(hub, candidate);
            }
        }

        // Assemblies defined by the hub record it in their genome json.  The hub is a managed hub only if it also
        // references genomes it does not define, and is not the current genome's own hub.
        Genome currentGenome = GenomeManager.getInstance().getCurrentGenome();
        boolean referencesHosted = genomeIds.stream().anyMatch(id -> !hub.hasTwoBit(id));
        if (referencesHosted && !isGenomeHub(currentGenome, hub)) {
            HubRegistry.addUserHub(hub);
        }

        if (currentGenome != null && hub.referencesGenome(currentGenome.getUCSCId())) {
            IGVMenuBar.getInstance().updateMenus(currentGenome);
            selectTracks(hub, currentGenome);
        } else {
            // Offer to switch to one of the hub's genomes now in the genome list.
            Map<String, GenomeListItem> updatedList = genomeListManager.getGenomeItemMap();
            List<GenomeListItem> hubGenomes = genomeIds.stream()
                    .map(id -> findInGenomeList(updatedList, id))
                    .filter(Objects::nonNull)
                    .toList();
            if (!hubGenomes.isEmpty()) {
                GenomeListItem choice =
                        HubGenomeSelectionDialog.selectGenomeToSwitchTo(IGV.getInstance().getMainFrame(), hubGenomes);
                if (choice != null) {
                    // No annotation prompt: the hub's track dialog is the one chooser
                    GenomeManager.getInstance().loadGenomeById(choice.getId(), false, true);
                    selectTracksForCurrentGenome(hub);
                }
            }
        }
    }

    /**
     * Create the local genome json for a genome selected from a multi-genome hub and add it to the genome list.
     * Default annotations are selected the first time the genome is loaded, so the json is not complete until then.
     */
    private static GenomeListItem addGenome(Hub hub, HubGenomeSelectionDialog.Candidate candidate) throws IOException {

        File localFile;
        if (candidate.config() != null) {
            // Assembly defined by the hub.  The hub becomes the genome's annotation hub (hubs[0]).
            GenomeConfig config = candidate.config();
            config.setHubs(Arrays.asList(hub.getUrl()));
            config.selectAnnotationsOnLoad = true;
            localFile = GenomeDownloadUtils.saveLocalGenome(config);
        } else {
            localFile = GenomeSelectionDialog.createLocalGenome(candidate.item(), false, false, false);
        }

        String path = localFile != null ? localFile.getAbsolutePath() : candidate.item().getPath();
        GenomeListItem item = new GenomeListItem(candidate.name(), path, candidate.id());
        GenomeListManager.getInstance().addGenomeItem(item);
        return item;
    }

    /**
     * True if the hub is one of the genome's own hubs (from its genome json), e.g. the hub the assembly was created
     * from.  Such a hub is listed with the genome's hubs and is not also a managed hub.
     */
    private static boolean isGenomeHub(Genome genome, Hub hub) {
        return genome != null && genome.getTrackHubs().stream().anyMatch(h -> h.getUrl().equals(hub.getUrl()));
    }

    /**
     * The genome list record for a hub genome id, which may name a UCSC hosted genome by its GenArk accession
     * (e.g. GCF_000001405.40 for hg38), or null if the genome is not in the list.
     */
    private static GenomeListItem findInGenomeList(Map<String, GenomeListItem> genomeList, String genomeId) {
        GenomeListItem item = genomeList.get(genomeId);
        if (item == null) {
            // The hub names a UCSC genome by database name or accession; a list entry may carry a GRC or 1kg name
            // for the same genome (GRCh38 for hg38), so compare by canonical UCSC name.
            String ucscId = Hub.getUcscIdForAlias(genomeId);
            String target = ucscId != null ? ucscId : genomeId;
            item = genomeList.get(target);
            if (item == null) {
                item = genomeList.values().stream()
                        .filter(i -> target.equals(Genome.toUCSCId(i.getId())))
                        .findFirst().orElse(null);
            }
        }
        return item;
    }

    /**
     * Return the genome list record for a genome id known to IGV -- already in the genome list, or hosted (including
     * GenArk) -- or null if the id is unknown.
     */
    private static GenomeListItem findKnownGenome(String genomeId) throws IOException {
        GenomeListItem item = findInGenomeList(GenomeListManager.getInstance().getGenomeItemMap(), genomeId);
        return item != null ? item : findHostedGenome(genomeId);
    }

    /**
     * The hosted genome record for a hub genome id.  IGV's own hosted list is keyed by UCSC name, so an accession is
     * tried under that name first; otherwise it would resolve to the GenArk record for the same assembly.
     */
    private static GenomeListItem findHostedGenome(String genomeId) {
        String ucscId = Hub.getUcscIdForAlias(genomeId);
        GenomeListItem item = ucscId == null ? null : HostedGenomes.getGenomeListItem(ucscId);
        return item != null ? item : HostedGenomes.getGenomeListItem(genomeId);
    }

    /**
     * Load a known genome as if it had been selected from the genome list, or from "Load Genome from Host" if it
     * is not in the list yet, without the default annotation prompt.
     */
    private static void loadKnownGenome(GenomeListItem item) throws IOException {
        String genome = item.getId();
        if (GenomeListManager.getInstance().getGenomeTableRecord(genome) == null) {
            File localFile = GenomeSelectionDialog.createLocalGenome(item, false, false, false);
            genome = localFile != null ? localFile.getAbsolutePath() : item.getPath();
        }
        GenomeManager.getInstance().loadGenomeById(genome, false, true);
    }

    private static void selectTracksForCurrentGenome(Hub hub) {
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        if (genome != null && hub.referencesGenome(genome.getUCSCId())) {
            selectTracks(hub, genome);
        }
    }

    private static void selectTracks(Hub hub, Genome genome) {
        String ucscId = genome.getUCSCId();
        if (hub.getSupportedTrackCount(ucscId) > 0) {
            SelectHubTracksAction.selectAndLoadTracks(hub, ucscId);
        } else {
            MessageUtils.showMessage(hub.getShortLabel() + " does not have any IGV supported tracks for " + ucscId);
        }
    }
}
