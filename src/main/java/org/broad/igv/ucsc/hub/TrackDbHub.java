package org.broad.igv.ucsc.hub;

import org.broad.igv.feature.genome.load.TrackConfig;
import org.broad.igv.ui.IGV;

import java.util.*;
import java.util.stream.Collectors;

import static org.broad.igv.ucsc.hub.Hub.getPriority;

public class TrackDbHub {

    List<Stanza> trackStanzas;
    List<Stanza> groupStanzas;

    static Set supportedTypes = new HashSet(Arrays.asList("bigBed", "bigWig", "bigGenePred", "vcfTabix", "refgene"));

    static Set filterTracks = new HashSet(Arrays.asList("cytoBandIdeo", "assembly", "gap", "gapOverlap", "allGaps",
            "cpgIslandExtUnmasked", "windowMasker"));

    static Map<String, String> vizModeMap = Map.of(
            "pack", "EXPANDED",
            "full", "EXPANDED",
            "squish", "SQUISHED",
            "dense", "COLLAPSED");


    public TrackDbHub(List<Stanza> trackStanzas, List<Stanza> groupStanzas) {

        this.groupStanzas = groupStanzas;
        this.trackStanzas = trackStanzas;
    }

    public List<TrackConfigGroup> getGroupedTrackConfigurations() {

        // Build map of group objects
        Map<String, TrackConfigGroup> groupMap = new HashMap<>();
        // Create a group for tracks without a group
        groupMap.put("", new TrackConfigGroup("", "", 0, true));

        if (this.groupStanzas != null) {
            for (Stanza groupStanza : this.groupStanzas) {
                String name = groupStanza.getProperty("name");
                boolean defaultOpen = "0".equals(groupStanza.getProperty("defaultIsClosed"));
                int priority = groupStanza.hasProperty("priority") ? getPriority(groupStanza.getProperty("priority")) : Integer.MAX_VALUE - 1;
                groupMap.put(name, new TrackConfigGroup(name, groupStanza.getProperty("label"), priority, defaultOpen));
            }
        }

        // Build map of stanzas to resolve parents
        Map<String, Stanza> trackStanzaMap = new HashMap<>();
        for (Stanza s : this.trackStanzas) {
            trackStanzaMap.put(s.getProperty("track"), s);
        }

        // Initialized cache of track containers.  Use linked hashmap to maintain insertion order
        LinkedHashMap<String, TrackConfigGroup> parentCache = new LinkedHashMap<>();

        // Tracks
        final List<TrackConfig> trackConfigs = this.getTrackConfigs(stanza -> !filterTracks.contains(stanza.name));
        for (TrackConfig c : trackConfigs) {
            String groupName = c.getGroup() != null ? c.getGroup() : "";

            // Some heads (at least one, the washu epigenomics hub) reference groups that are not defined.
            if (!groupMap.containsKey(groupName)) {
                groupMap.put(groupName, new TrackConfigGroup(groupName, groupName, Integer.MAX_VALUE, false));
            }

            final TrackConfigGroup trackConfigGroup = groupMap.get(groupName);
            int priority = trackConfigGroup.priority;
            trackConfigGroup.tracks.add(c);

            String parentName = c.getStanzaParent();
            if (parentName != null && trackStanzaMap.containsKey(parentName)) {
                // Create a contingent container, will be used if a sufficient # of tracks belong to this container
                TrackConfigGroup container = parentCache.get(parentName);
                if (container == null) {
                    Stanza s = trackStanzaMap.get(parentName);
                    String label = trackConfigGroup.label.equals("") ?
                            s.getProperty("shortLabel") :
                            trackConfigGroup.label + " - " + s.getProperty("shortLabel");
                    int p = s.hasProperty("priority") ? getPriority(s.getProperty("priority")) : priority + 1;
                    container = new TrackConfigGroup(parentName, label, p, false);
                    parentCache.put(parentName, container);
                }
                container.tracks.add(c);
            }
        }

        // Flatten the track groups into a list.  Remove empty groups.
        List<TrackConfigGroup> groupedTrackConfigurations = groupMap.values().stream()
                .filter(g -> !g.isEmpty()).collect(Collectors.toList());

        // Promote contingent embedded track groups (e.g. composite tracks) to top level group if # of tracks > threshold
        // Member tracks must also be removed from existing top level categories
        List<TrackConfigGroup> tmp = new ArrayList<>();
        Set<TrackConfig> toRemove = new HashSet<>();
        for (TrackConfigGroup parent : parentCache.values()) {
            if (parent.tracks.size() > 5) {
                tmp.add(parent);
                toRemove.addAll(parent.tracks);
            }
        }
        if (toRemove.size() > 0) {
            for (TrackConfigGroup g : groupedTrackConfigurations) {
                g.tracks = g.tracks.stream().filter(t -> !toRemove.contains(t)).collect(Collectors.toList());
            }
        }

        // Filter empty groups
        groupedTrackConfigurations = groupedTrackConfigurations.stream().filter(t -> t.tracks.size() > 0).collect(Collectors.toList());
        groupedTrackConfigurations.addAll(tmp);
        Collections.sort(groupedTrackConfigurations, Comparator.comparingInt(o -> o.priority));
        return groupedTrackConfigurations;
    }

    /**
     * Return an array of igv track config objects that satisfy the filter
     */
    List<TrackConfig> getTrackConfigs(java.util.function.Function<Stanza, Boolean> filter) {
        return this.trackStanzas.stream().filter(t -> {
                    return supportedTypes.contains(t.format()) && t.hasProperty("bigDataUrl") && (filter == null || filter.apply(t));
                })
                .map(t -> this.getTrackConfig(t))
                .toList();
    }


    private TrackConfig getTrackConfig(Stanza t) {

        String format = t.format();
        String url = t.getProperty("bigDataUrl");
        TrackConfig config = new TrackConfig(url);

        config.setPanelName(IGV.DATA_PANEL_NAME);

        config.setId(t.getProperty("track"));
        config.setName(t.getProperty("shortLabel"));

        // A rather complex workaround for some composite tracks
        String longLabel = t.getOwnProperty("longLabel");
        if(longLabel == null) {
            String inheritedLongLabel = t.getProperty("longLabel");
            longLabel = inheritedLongLabel.length() > config.getName().length() ? inheritedLongLabel : config.getName();
        }
        config.setLongLabel(longLabel);

        // TODO -- work on recognizing big* formats
        // config.format = t.format();

        config.setUrl(t.getProperty("bigDataUrl"));

        // Expanded display mode does not work well in IGV desktop for some tracks
        //config.displayMode = t.displayMode();

        if(t.hasProperty("bigDataIndex")) {
            config.setIndexURL(t.getProperty("bigDataIndex"));
        }

        if (t.hasProperty("longLabel") && t.hasProperty("html")) {
            config.setDescription("<a target=\"_blank\" href=\"" + (t.getProperty("html") + "\">" + t.getProperty("longLabel") + "</a>"));
        } else if (t.hasProperty("longLabel")) {
            config.setDescription(t.getProperty("longLabel"));
        }

        if (t.hasProperty("html")) {
            config.setHtml(t.getProperty("html"));
        }

        String vizProperty = t.getProperty("visibility");

        if (vizProperty != null && vizModeMap.containsKey(vizProperty)) {
            config.setDisplayMode(vizModeMap.get(vizProperty));
        }

        if(t.hasProperty("maxWindowToDraw")) {
            config.setVisibilityWindow(Integer.parseInt(t.getProperty("maxWindowToDraw")));
        }

        boolean visibility = t.hasProperty("compositeTrack") ?
                "on".equals(t.getProperty("compositeTrack")) :
                !("hide".equals(vizProperty));

        config.setVisible(visibility);

        if (t.hasProperty("autoScale")) {
            config.setAutoscale(t.getProperty("autoScale").toLowerCase().equals("on"));
        }
        if (t.hasProperty("maxHeightPixels")) {
            String[] tokens = t.getProperty("maxHeightPixels").split(":");
            config.setMaxHeight(Integer.parseInt(tokens[0]));
            config.setHeight(Integer.parseInt(tokens[1]));
            config.setMinHeight(Integer.parseInt(tokens[2]));
        }
        // TODO -- graphTypeDefault
        // TODO -- windowingFunction
        if (t.hasProperty("color")) {
            String c = t.getProperty("color");
            config.setColor(c.indexOf(",") > 0 ? "rgb(" + c + ")" : c);
        }
        if (t.hasProperty("altColor")) {
            String c = t.getProperty("altColor");
            config.setAltColor(c.indexOf(",") > 0 ? "rgb(" + c + ")" : c);
        }
        if (t.hasProperty("viewLimits")) {
            String[] tokens = t.getProperty("viewLimits").split(":");
            config.setMin(Float.parseFloat(tokens[0]));
            if (tokens.length > 1) {
                config.setMax(Float.parseFloat(tokens[1]));
            }
        }
        if (t.hasProperty("itemRgb")) {
            // TODO -- this not supported yet
        }
        if ("hide".equals(t.getProperty("visibility"))) {
            config.setVisible(false);
        }
        if (t.hasProperty("url")) {
            config.setInfoURL(t.getProperty("url"));
        }
        if (t.hasProperty("searchIndex")) {
            config.setSearchIndex(t.getProperty("searchIndex"));
        }
        if (t.hasProperty("searchTrix")) {
            config.setTrixURL(t.getProperty("searchTrix"));
        }

        if (t.hasProperty("group")) {
            config.setGroup(t.getProperty("group"));
        }

        if (t.parent != null) {
            config.setStanzaParent(t.parent.name);
        }

        return config;
    }

}
