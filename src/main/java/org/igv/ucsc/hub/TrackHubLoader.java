package org.igv.ucsc.hub;

import org.igv.feature.genome.Genome;
import org.igv.feature.genome.GenomeDownloadUtils;
import org.igv.feature.genome.GenomeManager;
import org.igv.feature.genome.HostedGenomes;
import org.igv.feature.genome.load.GenomeConfig;
import org.igv.feature.genome.load.HubGenomeLoader;
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
 * from "Load Genome from Host" and then showing the hub's track selection dialog.
 * <p>
 * Multiple genome hub: the user is offered the referenced genomes not already in the genome list, and those selected
 * are added, with default annotation selection deferred to their first load.  Then, if the current genome is one of
 * the hub's genomes its track selection dialog is shown, otherwise the user is offered a switch to one of them,
 * followed by the hub's track selection dialog if the chosen genome is a hosted one.
 * <p>
 * In every case the hub is recorded in the HubRegistry, so it appears in the "Track Hubs" menu for each genome it
 * references.  Long running; must not be called on the event thread.
 */
public class TrackHubLoader {

    private static Logger log = LogManager.getLogger(TrackHubLoader.class);

    public static void loadHub(String url) throws IOException {

        Hub hub = HubParser.loadHub(url);
        List<String> genomeIds = hub.getGenomeIds();
        if (genomeIds.isEmpty()) {
            throw new IOException("Hub does not reference any genomes: " + url);
        }
        if (genomeIds.size() == 1) {
            loadSingleGenomeHub(hub, genomeIds.get(0));
        } else {
            loadMultiGenomeHub(hub, genomeIds);
        }
    }

    private static void loadSingleGenomeHub(Hub hub, String genomeId) throws IOException {

        Genome currentGenome = GenomeManager.getInstance().getCurrentGenome();

        if (currentGenome != null && hub.referencesGenome(currentGenome.getUCSCId())) {
            HubRegistry.addUserHub(hub);
            IGVMenuBar.getInstance().updateMenus(currentGenome);
            selectTracks(hub, currentGenome);

        } else if (hub.hasTwoBit(genomeId)) {
            // The hub defines the assembly.  The hub becomes the genome's annotation hub (hubs[0]).
            GenomeConfig config = hub.getGenomeConfig(genomeId);
            config.setHubs(Arrays.asList(hub.getUrl()));
            HubGenomeLoader.selectAnnotationTracks(hub, config);
            File genomeFile = GenomeDownloadUtils.saveLocalGenome(config);
            GenomeManager.getInstance().loadGenome(genomeFile.getAbsolutePath());

        } else {
            // A UCSC hosted genome referenced by id.  Loading it may prompt for its default annotations (first load
            // only); the hub's own tracks are then offered so they are not overlooked.
            GenomeListItem item = findKnownGenome(genomeId);
            if (item == null) {
                MessageUtils.showMessage("The genome referenced by this hub (" + genomeId + ") is not known to IGV.");
                return;
            }
            HubRegistry.addUserHub(hub);
            loadKnownGenome(item);
            selectTracksForCurrentGenome(hub);
        }
    }

    private static void loadMultiGenomeHub(Hub hub, List<String> genomeIds) throws IOException {

        GenomeListManager genomeListManager = GenomeListManager.getInstance();
        Map<String, GenomeListItem> genomeList = genomeListManager.getGenomeItemMap();

        // Referenced genomes not already in the genome list, and ids IGV does not know.
        List<HubGenomeSelectionDialog.Candidate> candidates = new ArrayList<>();
        List<String> unknownIds = new ArrayList<>();
        for (String genomeId : genomeIds) {
            if (genomeList.containsKey(genomeId)) {
                continue;
            }
            if (hub.hasTwoBit(genomeId)) {
                GenomeConfig config = hub.getGenomeConfig(genomeId);
                candidates.add(new HubGenomeSelectionDialog.Candidate(genomeId, config.getName(), config, null));
            } else {
                GenomeListItem item = HostedGenomes.getGenomeListItem(genomeId);
                if (item == null) {
                    unknownIds.add(genomeId);
                } else {
                    candidates.add(new HubGenomeSelectionDialog.Candidate(genomeId, item.getDisplayableName(), null, item));
                }
            }
        }

        if (!candidates.isEmpty() || !unknownIds.isEmpty()) {
            HubGenomeSelectionDialog dlg =
                    new HubGenomeSelectionDialog(IGV.getInstance().getMainFrame(), hub, candidates, unknownIds);
            UIUtilities.invokeAndWaitOnEventThread(() -> dlg.setVisible(true));
            if (!dlg.isCanceled()) {
                for (HubGenomeSelectionDialog.Candidate candidate : dlg.getSelectedCandidates()) {
                    addGenome(hub, candidate);
                }
            }
        }

        HubRegistry.addUserHub(hub);

        Genome currentGenome = GenomeManager.getInstance().getCurrentGenome();
        if (currentGenome != null && hub.referencesGenome(currentGenome.getUCSCId())) {
            IGVMenuBar.getInstance().updateMenus(currentGenome);
            selectTracks(hub, currentGenome);
        } else {
            // Offer to switch to one of the hub's genomes now in the genome list.
            Map<String, GenomeListItem> updatedList = genomeListManager.getGenomeItemMap();
            List<GenomeListItem> hubGenomes = genomeIds.stream()
                    .map(updatedList::get)
                    .filter(Objects::nonNull)
                    .toList();
            if (!hubGenomes.isEmpty()) {
                GenomeListItem choice =
                        HubGenomeSelectionDialog.selectGenomeToSwitchTo(IGV.getInstance().getMainFrame(), hubGenomes);
                if (choice != null) {
                    GenomeManager.getInstance().loadGenomeById(choice.getId());
                    // For a hosted genome offer the hub's tracks after the switch.  An assembly created from the
                    // hub has the hub as its annotation hub, so its default annotations are the hub's tracks.
                    if (!hub.hasTwoBit(choice.getId())) {
                        selectTracksForCurrentGenome(hub);
                    }
                }
            }
        }
    }

    /**
     * Create the local genome json for a genome selected from a multi-genome hub and add it to the genome list.
     * Default annotations are selected the first time the genome is loaded, so the json is not complete until then.
     */
    private static void addGenome(Hub hub, HubGenomeSelectionDialog.Candidate candidate) throws IOException {

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
        GenomeListManager.getInstance().addGenomeItem(new GenomeListItem(candidate.name(), path, candidate.id()));
    }

    /**
     * Return the genome list record for a genome id known to IGV -- already in the genome list, or hosted (including
     * GenArk) -- or null if the id is unknown.
     */
    private static GenomeListItem findKnownGenome(String genomeId) {
        GenomeListItem item = GenomeListManager.getInstance().getGenomeTableRecord(genomeId);
        return item != null ? item : HostedGenomes.getGenomeListItem(genomeId);
    }

    /**
     * Load a known genome as if it had been selected from the genome list, or from "Load Genome from Host" if it
     * is not in the list yet.
     */
    private static void loadKnownGenome(GenomeListItem item) throws IOException {
        if (GenomeListManager.getInstance().getGenomeTableRecord(item.getId()) != null) {
            GenomeManager.getInstance().loadGenomeById(item.getId());
        } else {
            File localFile = GenomeSelectionDialog.createLocalGenome(item, false, false, true);
            GenomeManager.getInstance().loadGenome(localFile != null ? localFile.getAbsolutePath() : item.getPath());
        }
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
