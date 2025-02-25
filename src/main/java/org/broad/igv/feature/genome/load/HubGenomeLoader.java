package org.broad.igv.feature.genome.load;

import org.apache.commons.io.FileUtils;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ucsc.hub.Hub;
import org.broad.igv.ucsc.hub.HubParser;
import org.broad.igv.ucsc.hub.TrackConfigContainer;
import org.broad.igv.ui.IGV;
import org.broad.igv.ucsc.hub.TrackHubSelectionDialog;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Loads a "genome" from a UCSC track hub
 */
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

        Hub hub = HubParser.loadAssemblyHub(this.hubURL);

        GenomeConfig config = hub.getGenomeConfig();

        // Potentially override default tracks from hub with user selections

        // Check previous selections for this hub first
        // TODO -- Maintain track order?
        String key = "hub:" + this.hubURL;
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

            int count = 0;
            for (TrackConfigContainer g : groupedTrackConfigurations) {
                count += g.tracks.size();
            }

            // If the total # of tracks is >= 20 filter to "Gene" groups, usually a single group
            List<TrackConfigContainer> filteredGroups = count < 20 ?
                    groupedTrackConfigurations :
                    groupedTrackConfigurations.stream().filter(g -> g.label.startsWith("Gene")).collect(Collectors.toList());


            TrackHubSelectionDialog dlg = new TrackHubSelectionDialog(hub, filteredGroups, IGV.getInstance().getMainFrame());
            dlg.setVisible(true);

            if(!dlg.isCanceled()) {
                List<TrackConfig> selectedTracks = dlg.getSelectedConfigs();
                config.setTracks(selectedTracks);

                // Remember selections in user preferences
                List<String> names = selectedTracks.stream().map((trackConfig) -> trackConfig.getName()).toList();
                PreferencesManager.getPreferences().put(key, String.join(",", names));
            }
        }

        Genome genome = new Genome(config);
        genome.setGenomeHub(hub);



        String genomeJson = config.toJson();
        File dir = DirectoryManager.getGenomeCacheDirectory();
        FileUtils.write(new File(dir, config.getId() + ".json"), genomeJson);

        return genome;


    }
}
