package org.igv.feature.genome.load;

import org.igv.Globals;
import org.igv.feature.genome.Genome;
import org.igv.feature.genome.GenomeDownloadUtils;
import org.igv.feature.genome.GenomeManager;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.PreferencesManager;
import org.igv.ucsc.hub.Hub;
import org.igv.ucsc.hub.HubParser;
import org.igv.ucsc.hub.TrackConfigContainer;
import org.igv.ui.IGV;
import org.igv.ucsc.hub.TrackSelectionDialog;
import org.igv.ui.genome.GenomeListManager;
import org.igv.ui.util.MessageUtils;
import org.igv.ui.util.UIUtilities;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Loads a "genome" from a UCSC track hub.  Hub URLs entered by the user are handled by
 * {@link org.igv.ucsc.hub.TrackHubLoader}; this class remains as the GenomeLoader for hub URL genome paths, as
 * recorded in sessions and the genome list.
 */
public class HubGenomeLoader extends GenomeLoader {

    static Logger log = LogManager.getLogger(HubGenomeLoader.class);

    /**
     * Returns true if the string looks like a UCSC hub.txt URL.
     */
    public static boolean isHubURL(String obj) {
        return obj != null && obj.toLowerCase().endsWith("hub.txt");       // <= very crude, perhaps read the first line
    }

    /**
     * Converts a GenBank/RefSeq assembly accession (GCA_/GCF_) to its UCSC hub.txt URL, or null if not an accession.
     */
    public static String convertToHubURL(String accension) {
        //https://hgdownload.soe.ucsc.edu/hubs/GCF/016/808/095/GCF_016808095.1/
        //https://hgdownload.soe.ucsc.edu/hubs/GCA/028/534/965/GCA_028534965.1/
        if (accension.startsWith("GCF") || accension.startsWith("GCA") && accension.length() >= 13) {
            String prefix = accension.substring(0, 3);
            String n1 = accension.substring(4, 7);
            String n2 = accension.substring(7, 10);
            String n3 = accension.substring(10, 13);
            return "https://hgdownload.soe.ucsc.edu/hubs/" + prefix + "/" + n1 + "/" + n2 + "/" + n3 + "/" + accension + "/hub.txt";
        } else {
            return null;
        }
    }

    /**
     * Sets the genome's annotation tracks, either from a legacy preference entry for this hub or from
     * user selections in the track selection dialog (interactive mode only).
     */
    public static void selectAnnotationTracks(Hub hub, GenomeConfig config) throws IOException {

        // Potentially override default tracks from hub with user selections

        // Check previous selections for this hub first -- this is deprecated, kept for backward compatibility.
        String key = "hub:" + hub.getUrl();
        final List<TrackConfigContainer> groupedTrackConfigurations = hub.getGroupedTrackConfigurations(config.id);

        if (PreferencesManager.getPreferences().hasExplicitValue(key)) {
            Set<String> selectedTrackNames = new HashSet<>(Arrays.asList(PreferencesManager.getPreferences().get(key).split(",")));
            List<TrackConfig> selectedTracks = groupedTrackConfigurations.stream()
                    .flatMap(group -> group.tracks.stream())
                    .filter(trackConfig -> selectedTrackNames.contains(trackConfig.name))
                    .collect(Collectors.toList());
            config.setTracks(selectedTracks);
            PreferencesManager.getPreferences().remove(key);   // Remove this deprecated preference
        }

        // If running in interactive mode opend dialog to set tracks.
        else if (IGV.hasInstance() && !Globals.isBatch() && !Globals.isHeadless() && !Globals.isTesting()) {

            TrackSelectionDialog dlg;
            try {
                final TrackSelectionDialog[] holder = new TrackSelectionDialog[1];
                UIUtilities.invokeAndWaitOnEventThread(() -> {
                    holder[0] = TrackSelectionDialog.getTrackHubSelectionDialog(hub, config.getUcscID(), null, true,
                            GenomeManager.SELECT_ANNOTATIONS_MESSAGE);
                    holder[0].setVisible(true);
                });
                dlg = holder[0];
            } catch (Exception e) {
                log.error("Error opening or using TrackHubSelectionDialog: " + e.getMessage(), e);
                MessageUtils.showMessage("Error selecting annotation tracks from " + hub.getUrl() + ": " + e.getMessage());
                return;
            }

            if (dlg.isCanceled()) {
                return;
            }

            List<TrackConfig> selectedTracks = dlg.getSelectedConfigs();
            config.setTracks(selectedTracks);
        }
    }


    /*********************************************/
    // GenomeLoad interface follows -- provided for backward compatibility (IGV versions up to 2.19.2).
    // The current IGV version uses the static loadGenome(hubURL) method.

    String hubURL;

    /**
     * Legacy loader instance for the given hub URL.
     */
    public HubGenomeLoader(String hubURL) {
        this.hubURL = hubURL;
    }

    /**
     * Legacy load: parses the hub, applies any track selections, saves a local genome json, and returns the Genome.
     *
     * @return
     * @throws IOException
     */
    @Override
    public Genome loadGenome() throws IOException {

        Hub hub = HubParser.loadHub(hubURL);

        List<GenomeConfig> configs = hub.getGenomeConfigs();
        if (configs.isEmpty()) {
            throw new IOException("Hub does not define a genome sequence (no twoBitPath): " + hubURL);
        }
        GenomeConfig config = configs.get(0);
        selectAnnotationTracks(hub, config);

        // Search for list of tracks for this hub in the preferences.  This was used prior to version 2.19.2, current
        // versions of IGV store track information in a genome json file.
        String key = "hub:" + this.hubURL;
        if (PreferencesManager.getPreferences().hasExplicitValue(key)) {
            List<TrackConfig> selectedTracks = new ArrayList<>();
            Set<String> selectedTrackNames = new HashSet<>(Arrays.asList(PreferencesManager.getPreferences().get(key).split(",")));
            List<TrackConfigContainer> trackConfigGroups = hub.getGroupedTrackConfigurations(config.getUcscID());
            for (TrackConfigContainer group : trackConfigGroups) {
                List<TrackConfig> trackConfigs = group.findTracks(trackConfig -> selectedTrackNames.contains(trackConfig.name));
                for (TrackConfig trackConfig : trackConfigs) {
                    if (selectedTrackNames.contains(trackConfig.name)) {
                        selectedTracks.add(trackConfig);
                    }
                }
            }
            config.setTracks(selectedTracks);
        }

        // Save genome json for future loads and remove preferences and remote reference (Genome is now defined by
        // local genome json, preferences and remote reference not needed)
        GenomeDownloadUtils.saveLocalGenome(config);
        GenomeListManager.getInstance().removeRemoteItem(config.id);

        return new Genome(config);
    }

}
