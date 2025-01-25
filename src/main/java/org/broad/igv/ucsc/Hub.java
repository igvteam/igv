package org.broad.igv.ucsc;

import org.broad.igv.Globals;
import org.broad.igv.feature.genome.load.GenomeConfig;
import org.broad.igv.feature.genome.load.TrackConfig;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;
import java.util.function.Function;

public class Hub {

    private final String url;
    String baseURL;
    Stanza hub;
    Stanza genomeStanza;
    List<Stanza> trackStanzas;
    List<Stanza> groupStanzas;
    Map<String, Integer> groupPriorityMap;

    static Set supportedTypes = new HashSet(Arrays.asList("bigBed", "bigWig", "bigGenePred", "vcfTabix"));
    static Set filterTracks = new HashSet(Arrays.asList("cytoBandIdeo", "assembly", "gap", "gapOverlap", "allGaps",
            "cpgIslandExtUnmasked", "windowMasker"));

    public static Hub loadHub(String url) throws IOException {

        int idx = url.lastIndexOf("/");
        String baseURL = url.substring(0, idx + 1);

        // Load stanzas
        List<Stanza> stanzas = loadStanzas(url);

        // Validation checks
        if (stanzas.size() < 2) {
            throw new RuntimeException("Expected at least 2 stanzas, hub and genome");
        }
        if (!"hub".equals(stanzas.get(0).type)) {
            throw new RuntimeException("Unexpected hub.txt file -- does the first line start with 'hub'?");
        }
        if (!"on".equals(stanzas.get(0).getProperty("useOneFile"))) {
            throw new RuntimeException("Only 'useOneFile' hubs are currently supported");
        }
        if (!"genome".equals(stanzas.get(1).type)) {
            throw new RuntimeException("Unexpected hub file -- expected 'genome' stanza but found " + stanzas.get(1).type);
        }

        // Load groups
        List<Stanza> groups = null;
        Stanza genomeStanza = stanzas.get(1);
        if (genomeStanza.hasProperty("groups")) {
            String groupsTxtURL = baseURL + genomeStanza.getProperty("groups");
            groups = loadStanzas(groupsTxtURL);
        }

        // load includes.  Nested includes are not supported
        List<Stanza> includes = stanzas.stream().filter(s -> "include".equals(s.type)).toList();
        for (Stanza s : includes) {
            List<Stanza> includeStanzas = loadStanzas(baseURL + s.getProperty("include"));
            for (Stanza inc : includeStanzas) {
                s.setProperty("visibility", "hide");
                stanzas.add(s);
            }
        }

        return new Hub(url, stanzas, groups);
    }

    private Hub(String url, List<Stanza> stanzas, List<Stanza> groupStanzas) {

        this.url = url;

        int idx = url.lastIndexOf("/");
        String baseURL = url.substring(0, idx + 1);
        this.baseURL = baseURL;

        if (stanzas.size() < 2) {
            throw new RuntimeException("Expected at least 2 stanzas, hub and genome");
        }
        // The first stanza must be type = hub
        if ("hub".equals(stanzas.get(0).type)) {
            this.hub = stanzas.get(0);
        } else {
            throw new RuntimeException("Unexpected hub.txt file -- does the first line start with 'hub'?");
        }
        if (!"on".equals(this.hub.getProperty("useOneFile"))) {
            throw new RuntimeException("Only 'useOneFile' hubs are currently supported");
        }


        // The second stanza should be a genome
        if ("genome".equals(stanzas.get(1).type)) {
            this.genomeStanza = stanzas.get(1);
        } else {
            throw new RuntimeException("Unexpected hub file -- expected 'genome' stanza but found " + stanzas.get(1).type);
        }

        // Remaining stanzas should be tracks
        this.trackStanzas = new ArrayList<>();
        for (int i = 2; i < stanzas.size(); i++) {
            if ("track".equals(stanzas.get(i).type)) {
                this.trackStanzas.add(stanzas.get(i));
            }
        }

        if (groupStanzas != null) {
            this.groupStanzas = groupStanzas;
            this.groupPriorityMap = new HashMap<>();
            for (Stanza g : groupStanzas) {
                if (g.hasProperty("priority")) {
                    this.groupPriorityMap.put(g.getProperty("name"), getPriority(g));
                }
            }
        }
    }

    /**
     * Retun the priority for the group.  The priority format is uncertain, but extends to at least 2 levels (e.g. 3.4).
     * Ignore levels > 2
     *
     * @param g the group stanza
     * @return A priority as an integer
     */
    private static int getPriority(Stanza g) {
        String priorityString = g.getProperty("priority");
        String [] tokens = priorityString.split("\\.");
        int p = Integer.parseInt(tokens[0]) * 10;
        if(tokens.length > 1) {
            p += Integer.parseInt(tokens[1]);
        }
        return p;
    }

    private static List<Hub.Stanza> loadStanzas(String url) throws IOException {
        List<Stanza> nodes = new ArrayList<>();
        Stanza currentNode = null;
        boolean startNewNode = true;
        int order = 0;
        try (BufferedReader br = ParsingUtils.openBufferedReader(url)) {
            String line;
            while ((line = br.readLine()) != null) {

                int indent = indentLevel(line);
                int i = line.indexOf(" ", indent);
                if (i < 0) {
                    // Break - start a new node
                    startNewNode = true;
                } else {
                    String key = line.substring(indent, i);
                    if (key.startsWith("#")) continue;
                    String value = line.substring(i + 1);
                    if (startNewNode) {
                        // Start a new node -- indent is currently ignored as igv.js does not support sub-tracks,
                        // so track stanzas are flattened
                        Stanza newNode = new Stanza(++order, key, value);
                        nodes.add(newNode);
                        currentNode = newNode;
                        startNewNode = false;
                    }
                    currentNode.setProperty(key, value);
                }
            }
        }
        return resolveParents(nodes);

    }

    private static int indentLevel(String str) {
        int level;
        for (level = 0; level < str.length(); level++) {
            char c = str.charAt(level);
            if (c != ' ' && c != '\t') break;
        }
        return level;
    }

    private static List<Stanza> resolveParents(List<Stanza> nodes) {
        Map<String, Stanza> nodeMap = new HashMap<>();
        for (Stanza n : nodes) {
            nodeMap.put(n.name, n);
        }
        for (Stanza n : nodes) {
            if (n.properties.containsKey("parent")) {
                String parentName = firstWord(n.properties.get("parent"));
                n.parent = nodeMap.get(parentName);
            }
        }
        return nodes;
    }

    public GenomeConfig getGenomeConfig(boolean includeTracks) {
        // TODO -- add blat?  htmlPath?

        GenomeConfig config = new GenomeConfig();
        config.setId(this.genomeStanza.getProperty("genome"));
        if (this.genomeStanza.hasProperty("scientificName")) {
            config.setName(this.genomeStanza.getProperty("scientificName"));
        } else if (this.genomeStanza.hasProperty("organism")) {
            config.setName(this.genomeStanza.getProperty("organism"));
        } else if (this.genomeStanza.hasProperty("description")) {
            config.setName(this.genomeStanza.getProperty("description"));
        }
        if (config.getName() == null) {
            config.setName(config.getId());
        } else {
            config.setName(config.getName() + " (" + config.getId() + ")");
        }

        config.setTwoBitURL(this.baseURL + this.genomeStanza.getProperty("twoBitPath"));
        config.setNameSet("ucsc");
        config.setWholeGenomeView(false);

        if (this.genomeStanza.hasProperty("defaultPos")) {
            config.setDefaultPos(this.genomeStanza.getProperty("defaultPos"));
        }

        config.setDescription(config.getId());

        if (this.genomeStanza.hasProperty("blat")) {
            config.setBlat(this.baseURL + this.genomeStanza.getProperty("blat"));
        }
        if (this.genomeStanza.hasProperty("chromAliasBb")) {
            config.setChromAliasBbURL(this.baseURL + this.genomeStanza.getProperty("chromAliasBb"));
        }
        if (this.genomeStanza.hasProperty("twoBitBptURL")) {
            config.setTwoBitBptURL(this.baseURL + this.genomeStanza.getProperty("twoBitBptURL"));
        }

        if (this.genomeStanza.hasProperty("twoBitBptUrl")) {
            config.setTwoBitBptURL(this.baseURL + this.genomeStanza.getProperty("twoBitBptUrl"));
        }

        // chromSizes can take a very long time to load, and is not useful with the default WGV = off
        // if (this.genomeStanza.hasProperty("chromSizes")) {
        //     config.chromSizes = this.baseURL + this.genomeStanza.getProperty("chromSizes")
        // }

        if (this.genomeStanza.hasProperty("description")) {
            config.setDescription(config.getDescription() + "\n" + this.genomeStanza.getProperty("description"));
        }
        if (this.genomeStanza.hasProperty("organism")) {
            config.setDescription(config.getDescription() + "\n" + this.genomeStanza.getProperty("organism"));
        }
        if (this.genomeStanza.hasProperty("scientificName")) {
            config.setDescription(config.getDescription() + "\n" + this.genomeStanza.getProperty("scientificName"));
        }

        if (this.genomeStanza.hasProperty("htmlPath")) {
            config.setInfoURL(this.baseURL + this.genomeStanza.getProperty("htmlPath"));
        }

        // Search for cytoband
        /*
        track cytoBandIdeo
        shortLabel Chromosome Band (Ideogram)
        longLabel Ideogram for Orientation
        group map
        visibility dense
        type bigBed 4 +
        bigDataUrl bbi/GCA_004027145.1_DauMad_v1_BIUU.cytoBand.bb
         */
        for (Stanza t : this.trackStanzas) {
            if ("cytoBandIdeo".equals(t.name) && t.hasProperty("bigDataUrl")) {
                config.setCytobandBbURL(this.baseURL + t.getProperty("bigDataUrl"));
                break;
            }
        }

        // Tracks.  To prevent loading tracks set `includeTrackGroups`to false or "none"
        if (includeTracks) {
            Function<Stanza, Boolean> filter = (t) -> {
                return !Hub.filterTracks.contains(t.name) &&
                        (!"hide".equals(t.getProperty("visibility")));
            };
            config.setTracks(this.getTracksConfigs(filter));
        }

        // config.trackConfigurations = this.#getGroupedTrackConfigurations()

        return config;
    }

    public List<TrackConfigGroup> getGroupedTrackConfigurations() {

        // Organize track configs by group
        LinkedHashMap<String, List<TrackConfig>> trackConfigMap = new LinkedHashMap<>();
        java.util.function.Function<Stanza, Boolean> filter = (stanza -> !stanza.name.equals("cytoBandIdeo"));
        for (TrackConfig c : this.getTracksConfigs(filter)) {
            String groupName = c.getGroup() != null ? c.getGroup() : "other";
            if (!trackConfigMap.containsKey(groupName)) {
                trackConfigMap.put(groupName, new ArrayList<>());
            }
            trackConfigMap.get(groupName).add(c);
        }

        // Extract map of group names
        Map<String, String> groupNamesMap = new HashMap<>();
        if (this.groupStanzas != null) {
            for (Stanza groupStanza : this.groupStanzas) {
                groupNamesMap.put(groupStanza.getProperty("name"), groupStanza.getProperty("label"));
            }
        }

        // Use linked has map to maintain order
        List<TrackConfigGroup> groupedTrackConfigurations = new ArrayList<>();
        for (Map.Entry<String, List<TrackConfig>> entry : trackConfigMap.entrySet()) {
            String group = entry.getKey();
            String label = groupNamesMap.containsKey(group) ? groupNamesMap.get(group) : group;
            groupedTrackConfigurations.add(new TrackConfigGroup(label, entry.getValue()));

        }
        return groupedTrackConfigurations;
    }

    /**
     * Return an array of igv track config objects that satisfy the filter
     */
    List<TrackConfig> getTracksConfigs(java.util.function.Function<Stanza, Boolean> filter) {
        return this.trackStanzas.stream().filter(t -> {
                    return supportedTypes.contains(t.format()) && t.hasProperty("bigDataUrl") && (filter == null || filter.apply(t));
                })
                .map(t -> this.getTrackConfig(t))
                .toList();
    }

    TrackConfig getTrackConfig(Stanza t) {

        String format = t.format();
        String url = this.baseURL + t.getProperty("bigDataUrl");
        TrackConfig config = new TrackConfig(url);

        config.setPanelName(IGV.DATA_PANEL_NAME);

        config.setId(t.getProperty("track"));
        config.setName(t.getProperty("shortLabel"));

        // TODO -- work on recognizing big* formats
        // config.format = t.format();

        config.setUrl(this.baseURL + t.getProperty("bigDataUrl"));

        // Expanded display mode does not work well in IGV desktop for some tracks
        //config.displayMode = t.displayMode();

        if ("vcfTabix".equals(format)) {
            config.setIndexURL(config.getUrl() + ".tbi");
        }

        if (t.hasProperty("longLabel") && t.hasProperty("html")) {
            config.setDescription("<a target=\"_blank\" href=\"" + (this.baseURL + t.getProperty("html")) + "\">" + t.getProperty("longLabel") + "</a>");
        } else if (t.hasProperty("longLabel")) {
            config.setDescription(t.getProperty("longLabel"));
        }

        if (t.hasProperty("html")) {
            config.setHtml(this.baseURL + t.getProperty("html"));
        }

        config.setVisible(!("hide".equals(t.getProperty("visibility"))));

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
            config.setTrixURL(this.baseURL + t.getProperty("searchTrix"));
        }

        if (t.hasProperty("group")) {
            config.setGroup(t.getProperty("group"));
            if (this.groupPriorityMap != null && this.groupPriorityMap.containsKey(config.getGroup())) {
                int nextPriority = this.groupPriorityMap.get(config.getGroup()) + 1;
                config.setOrder(nextPriority);
                this.groupPriorityMap.put(config.getGroup(), nextPriority);
            }
        }

        return config;
    }

    public String getUrl() {
        return url;
    }


    static class Stanza {

        private final String type;
        private final String name;
        private final int order;

        private Map<String, String> properties;

        Stanza parent;

        Stanza(int order, String type, String name) {
            this.order = order;
            this.type = type;
            this.name = name;
            this.properties = new HashMap<>();
        }

        void setProperty(String key, String value) {
            this.properties.put(key, value);
        }

        String getProperty(String key) {
            if (this.properties.containsKey(key)) {
                return this.properties.get(key);
            } else if (this.parent != null) {
                return this.parent.getProperty(key);
            } else {
                return null;
            }
        }

        boolean hasProperty(String key) {
            if (this.properties.containsKey(key)) {
                return true;
            } else if (this.parent != null) {
                return this.parent.hasProperty(key);
            } else {
                return false;
            }
        }

        String format() {
            String type = this.getProperty("type");
            if (type != null) {
                // Trim extra bed qualifiers (e.g. bigBed + 4)
                return firstWord(type);
            }
            return null;  // unknown type
        }

        /**
         * IGV display mode
         */
        String displayMode() {
            String viz = this.getProperty("visibility");
            if (viz == null) {
                return "COLLAPSED";
            } else {
                viz = viz.toLowerCase();
                switch (viz) {
                    case "dense":
                        return "COLLAPSED";
                    case "pack":
                        return "EXPANDED";
                    case "squish":
                        return "SQUISHED";
                    default:
                        return "COLLAPSED";
                }
            }
        }
    }

    static String firstWord(String str) {
        return Globals.whitespacePattern.split(str)[0];
    }

}
