package org.broad.igv.ucsc.hub;

import org.broad.igv.feature.genome.load.TrackConfig;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.function.Function;

public class TrackConfigContainer {

    public int priority;
    public String name;
    public String label;
    public boolean defaultOpen;
    public List<TrackConfig> tracks;
    public List<TrackConfigContainer> children;

    public TrackConfigContainer(String name, String label, int priority, boolean defaultOpen) {
        this.name = name;
        this.priority = priority;
        this.label = label;
        this.defaultOpen = defaultOpen;
        this.tracks = new ArrayList<>();
        this.children = new ArrayList<>();
    }

    public boolean isEmpty() {
        return tracks.isEmpty() && (children == null || children.isEmpty() || children.stream().allMatch(TrackConfigContainer::isEmpty));
    }

    public void map(Function<TrackConfig, Void> f) {
        for (TrackConfig config : tracks) {
            f.apply(config);
        }
        for (TrackConfigContainer g : children) {
            g.map(f);
        }
    }

    public void findSelectedConfigs(List<TrackConfig> selectedConfigs) {
        for (TrackConfig trackConfig : selectedConfigs) {
            if (trackConfig.getVisible() == true) {
                selectedConfigs.add(trackConfig);
            }
        }
        for (TrackConfigContainer container : children) {
            container.findSelectedConfigs(selectedConfigs);
        }
    }

    public int countTracks() {
        int count = tracks.size();
        for (TrackConfigContainer container : children) {
            count += container.countTracks();
        }
        return count;
    }

    public int countSelectedTracks() {
        int count = 0;
        for (TrackConfig trackConfig : tracks) {
            if (trackConfig.getVisible() == true) {
                count++;
            }
        }
        for (TrackConfigContainer container : children) {
            count += container.countSelectedTracks();
        }
        return count;
    }

    public void trim() {
        children.stream().filter(c -> !c.isEmpty()).forEach(c -> c.trim());
    }

    public void setTrackVisibility(Set<String> loadedTrackPaths) {
        for (TrackConfig trackConfig : tracks) {
            trackConfig.setVisible(loadedTrackPaths.contains(trackConfig.getUrl()));
        }
        for (TrackConfigContainer container : children) {
            container.setTrackVisibility(loadedTrackPaths);
        }
    }
}
