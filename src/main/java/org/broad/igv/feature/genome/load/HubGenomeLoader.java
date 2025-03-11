package org.broad.igv.feature.genome.load;

import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeDownloadUtils;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ucsc.hub.Hub;
import org.broad.igv.ucsc.hub.HubParser;
import org.broad.igv.ucsc.hub.TrackConfigContainer;
import org.broad.igv.ui.IGV;
import org.broad.igv.ucsc.hub.TrackHubSelectionDialog;
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.ui.commandbar.GenomeListManager;
import org.broad.igv.ui.util.MessageUtils;

import javax.swing.*;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Loads a "genome" from a UCSC track hub
 */
public class HubGenomeLoader extends GenomeLoader {

    static Logger log = LogManager.getLogger(HubGenomeLoader.class);

    public static boolean isHubURL(String obj) {
        return obj != null && obj.toLowerCase().endsWith("hub.txt");       // <= very crude, perhaps read the first line
    }

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
     * Load and parse a hub.txt file.  Creates and stores a genome json definition for future loads.
     *
     * @param hubURL
     * @return
     * @throws IOException
     */
    public static void loadGenome(String hubURL) throws IOException {

        final WaitCursorManager.CursorToken token = WaitCursorManager.showWaitCursor();

        SwingWorker<File, Void> worker = new SwingWorker<>() {
            @Override
            protected File doInBackground() throws Exception {
                File genomeFile = null;
                try {
                    Hub hub = HubParser.loadAssemblyHub(hubURL);
                    final GenomeConfig config = getGenomeConfig(hub);
                    genomeFile = GenomeDownloadUtils.saveLocalGenome(config);
                } catch (Exception e) {
                    log.error("Error loading hub: " + e.getMessage());
                    MessageUtils.showMessage("Error loading hub: " + e.getMessage());
                }
                return genomeFile;
            }

            @Override
            protected void done() {
                try {
                    File genomeFile = get();
                    WaitCursorManager.removeWaitCursor(token);
                    GenomeManager.getInstance().loadGenome(genomeFile.getAbsolutePath());
                } catch (Exception e) {
                    log.error("Error loading hub: " + e.getMessage());
                    MessageUtils.showMessage("Error loading hub: " + e.getMessage());
                } finally {
                    WaitCursorManager.removeWaitCursor(token);
                }
            }
        };

        worker.execute();
    }

    /**
     * Shortcut method provided to support load-from-url.  The hub must be read and parsed to determine if it is an
     * assembly hub.  This method avoids the need to load and parse again.
     *
     * @param hub
     * @throws IOException
     */
    public static void loadAssemblyHub(Hub hub) throws IOException {
        final GenomeConfig config = getGenomeConfig(hub);
        File genomeFile = GenomeDownloadUtils.saveLocalGenome(config);
        GenomeManager.getInstance().loadGenome(genomeFile.getAbsolutePath());
    }


    private static GenomeConfig getGenomeConfig(Hub hub) throws IOException {

        GenomeConfig config = hub.getGenomeConfig();

        // Potentially override default tracks from hub with user selections

        // Check previous selections for this hub first
        // TODO -- Maintain track order?
        String key = "hub:" + hub.getUrl();
        final List<TrackConfigContainer> groupedTrackConfigurations = hub.getGroupedTrackConfigurations();

        if (PreferencesManager.getPreferences().hasExplicitValue(key)) {
            Set<String> selectedTrackNames = new HashSet<>(Arrays.asList(PreferencesManager.getPreferences().get(key).split(",")));
            List<TrackConfig> selectedTracks = groupedTrackConfigurations.stream()
                    .flatMap(group -> group.tracks.stream())
                    .filter(trackConfig -> selectedTrackNames.contains(trackConfig.getName()))
                    .collect(Collectors.toList());
            config.setTracks(selectedTracks);
        }

        // If running in interactive mode opend dialog to set tracks.
        else if (IGV.hasInstance() && !Globals.isBatch() && !Globals.isHeadless() && !Globals.isTesting()) {

            TrackHubSelectionDialog dlg = TrackHubSelectionDialog.getTrackHubSelectionDialog(hub, null, true);

            boolean dlgSuccess = true;
            try {
                SwingUtilities.invokeAndWait(() -> dlg.setVisible(true));
            } catch (Exception e) {
                dlgSuccess = false;
                log.error("Error opening or using TrackHubSelectionDialog: " + e.getMessage());
            }


            if (!dlg.isCanceled() && dlgSuccess) {
                List<TrackConfig> selectedTracks = dlg.getSelectedConfigs();
                config.setTracks(selectedTracks);

                // Remember selections in user preferences
                // List<String> names = selectedTracks.stream().map((trackConfig) -> trackConfig.getName()).toList();
                // PreferencesManager.getPreferences().put(key, String.join(",", names));
                // TODO -- read these for backward compatibility?
            }
        }
        config.setHubs(Arrays.asList(hub.getUrl()));
        Genome genome = new Genome(config);
        genome.setGenomeHub(hub);
        return config;
    }


    /*********************************************/
    // GenomeLoad interface follows -- provided for backward compatibility (IGV versions up to 2.19.2).
    // The current IGV version uses the static loadGenome(hubURL) method.

    String hubURL;

    public HubGenomeLoader(String hubURL) {
        this.hubURL = hubURL;
    }

    /**
     *
     * @return
     * @throws IOException
     */
    @Override
    public Genome loadGenome() throws IOException {

        Hub hub = HubParser.loadAssemblyHub(hubURL);

        GenomeConfig config = getGenomeConfig(hub);

        // Search for list of tracks for this hub in the preferences.  This was used prior to version 2.19.2, current
        // versions of IGV store track information in a genome json file.
        String key = "hub:" + this.hubURL;
        if (PreferencesManager.getPreferences().hasExplicitValue(key)) {
            List<TrackConfig> selectedTracks = new ArrayList<>();
            Set<String> selectedTrackNames = new HashSet<>(Arrays.asList(PreferencesManager.getPreferences().get(key).split(",")));
            List<TrackConfigContainer> trackConfigGroups = hub.getGroupedTrackConfigurations();
            for (TrackConfigContainer group : trackConfigGroups) {
                List<TrackConfig> trackConfigs = group.findTracks(trackConfig -> selectedTrackNames.contains(trackConfig.getName()));
                for (TrackConfig trackConfig : trackConfigs) {
                    if (selectedTrackNames.contains(trackConfig.getName())) {
                        selectedTracks.add(trackConfig);
                    }
                }
            }
            config.setTracks(selectedTracks);
        }

        // Save genome json for future loads and remove preferences and remote reference (Genome is now defined by
        // local genome json, preferences and remote reference not needed)
        GenomeDownloadUtils.saveLocalGenome(config);
        PreferencesManager.getPreferences().remove(key);
        GenomeListManager.getInstance().removeRemoteItem(config.getId());


        Genome genome = new Genome(config);
        genome.setGenomeHub(hub);
        return genome;
    }

}
