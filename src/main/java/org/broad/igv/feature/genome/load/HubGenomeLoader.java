package org.broad.igv.feature.genome.load;

import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ucsc.Hub;
import org.broad.igv.ucsc.TrackConfigGroup;
import org.broad.igv.ui.IGV;
import org.broad.igv.ucsc.HubTrackSelectionDialog;

import java.io.IOException;
import java.util.*;

public class HubGenomeLoader extends GenomeLoader {

    String hubURL;

    public HubGenomeLoader(String hubURL) {
        this.hubURL = hubURL;
    }

    public static boolean isHubURL(String obj) {
        return obj.endsWith("/hub.txt");
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

    @Override
    public Genome loadGenome() throws IOException {

        Hub hub = Hub.loadHub(this.hubURL);

        GenomeConfig config = hub.getGenomeConfig(false);

        // Potentially override default tracks from hub with user selections

        // Check previous selections for this hub first
        // TODO -- Maintain track order?
        String key = "hub:" + this.hubURL;
        if (PreferencesManager.getPreferences().hasExplicitValue(key)) {
            List<TrackConfig> selectedTracks = new ArrayList<>();
            Set<String> selectedTrackNames = new HashSet<>(Arrays.asList(PreferencesManager.getPreferences().get(key).split(",")));
            List<TrackConfigGroup> trackConfigGroups = hub.getGroupedTrackConfigurations();
            for (TrackConfigGroup group : trackConfigGroups) {
                for (TrackConfig trackConfig : group.tracks) {
                    if (selectedTrackNames.contains(trackConfig.getName())) {
                        selectedTracks.add(trackConfig);
                    }
                }
            }
            selectedTracks.sort((o1, o2) -> o1.getOrder() - o2.getOrder());
            config.setTracks(selectedTracks);
        }

        // If running in interactive mode opend dialog to set tracks.
        else if (IGV.hasInstance() && !Globals.isBatch() && !Globals.isHeadless() && !Globals.isTesting()) {
            HubTrackSelectionDialog dlg = new HubTrackSelectionDialog(hub.getGroupedTrackConfigurations(), IGV.getInstance().getMainFrame());
            dlg.setVisible(true);

            List<TrackConfig> selectedTracks = dlg.getSelectedConfigs();
            config.setTracks(selectedTracks);

            // Remember selections in user preferences
            List<String> names = selectedTracks.stream().map((trackConfig) -> trackConfig.getName()).toList();
            PreferencesManager.getPreferences().put(key, String.join(",", names));
        }

        Genome genome = new Genome(config);
        genome.setHub(hub);

        return genome;


    }
}
