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

    public List<TrackConfig> findTracks(Function<TrackConfig, Boolean> filter) {
        List<TrackConfig> found = new ArrayList<>();
        find(found, filter);
        return found;
    }

    private void find(List<TrackConfig> found, Function<TrackConfig, Boolean> filter) {
        for (TrackConfig config : tracks) {
            if (filter.apply(config)) {
                found.add(config);
            }
        }
        for (TrackConfigContainer g : children) {
            g.find(found, filter);
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
