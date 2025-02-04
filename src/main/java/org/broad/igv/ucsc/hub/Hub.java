package org.broad.igv.ucsc.hub;

import org.broad.igv.Globals;
import org.broad.igv.feature.genome.load.GenomeConfig;
import org.broad.igv.feature.genome.load.TrackConfig;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

public class Hub {

    private static Logger log = LogManager.getLogger(Hub.class);
    private final String url;
    private final String host;
    String baseURL;
    Stanza hubStanza;
    Stanza genomeStanza;
    List<Stanza> trackStanzas;
    List<Stanza> groupStanzas;
    Map<String, Integer> groupPriorityMap;

    static Set supportedTypes = new HashSet(Arrays.asList("bigBed", "bigWig", "bigGenePred", "vcfTabix"));
    static Set filterTracks = new HashSet(Arrays.asList("cytoBandIdeo", "assembly", "gap", "gapOverlap", "allGaps",
            "cpgIslandExtUnmasked", "windowMasker"));
    static Map<String, String> vizModeMap = Map.of(
            "pack", "EXPANDED",
            "full", "EXPANDED",
            "squish", "SQUISHED",
            "dense", "COLLAPSED");

    public static Hub loadHub(String url) throws IOException {
        return new Hub(url);
    }

    Hub(String url) throws IOException {

        this.url = url;

        int idx = url.lastIndexOf("/");
        String baseURL = url.substring(0, idx + 1);
        this.baseURL = baseURL;

        if (url.startsWith("https://") || url.startsWith("http://")) {
            try {
                URL tmp = new URL(url);
                this.host = tmp.getProtocol() + "://" + tmp.getHost();
            } catch (MalformedURLException e) {
                // This should never happen
                log.error("Error parsing base URL host", e);
                throw new RuntimeException(e);
            }
        } else {
            // Local file, no host
            this.host = "";
        }

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

        this.hubStanza = stanzas.get(0);

        // Load groups
        this.genomeStanza = stanzas.get(1);
        if (genomeStanza.hasProperty("groups")) {
            String groupsTxtURL = getDataURL(genomeStanza.getProperty("groups"));
            this.groupStanzas = loadStanzas(groupsTxtURL);
        }

        // load includes.  Nested includes (includes within includes) are not supported
        List<Stanza> includes = stanzas.stream().filter(s -> "include".equals(s.type)).toList();
        for (Stanza s : includes) {
            List<Stanza> includeStanzas = loadStanzas(getDataURL(s.getProperty("include")));
            for (Stanza inc : includeStanzas) {
                inc.setProperty("visibility", "hide");
                stanzas.add(inc);
            }
        }


        // The first stanza must be type = hub
        // The second stanza should be a genome
        // Remaining stanzas should be tracks
        this.trackStanzas = new ArrayList<>();
        for (int i = 2; i < stanzas.size(); i++) {
            if ("track".equals(stanzas.get(i).type)) {
                this.trackStanzas.add(stanzas.get(i));
            }
        }

        if (groupStanzas != null) {
            this.groupPriorityMap = new HashMap<>();
            for (Stanza g : groupStanzas) {
                if (g.hasProperty("priority")) {
                    this.groupPriorityMap.put(g.getProperty("name"), getPriority(g.getProperty("priority")));
                }
            }
        }
    }


    public String getShortLabel() {
        return this.hubStanza.hasProperty("shortLabel") ? this.hubStanza.getProperty("shortLabel") : this.url;
    }

    public String getLongLabel() {
        return this.hubStanza.hasProperty("longLabel") ? this.hubStanza.getProperty("longLabel") : this.url;
    }

    public String getDescriptionURL() {
        return this.hubStanza.hasProperty("descriptionUrl") ?
                this.getDataURL(this.hubStanza.getProperty("descriptionUrl")) :
                this.hubStanza.hasProperty("desriptionUrl") ?
                        this.getDataURL(this.hubStanza.getProperty("desriptionUrl")) : null;
    }


    /**
     * Return the priority for the group.  The priority format is uncertain, but extends to at least 2 levels (e.g. 3.4).
     * Ignore levels > 3
     *
     * @param priorityString Priority as a string (e.g. 3.4)
     * @return A priority as an integer
     */
    private static int getPriority(String priorityString) {
        String[] tokens = priorityString.split("\\.");
        int p = Integer.parseInt(tokens[0]) * 100;
        if (tokens.length > 1) {
            p += Integer.parseInt(tokens[1]) * 10;
        }
        if (tokens.length > 2) {
            p += Integer.parseInt(tokens[2]);
        }
        return p;
    }

    static List<Hub.Stanza> loadStanzas(String url) throws IOException {
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
                        Stanza newNode = new Stanza(key, value);
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

        config.setTwoBitURL(getDataURL(this.genomeStanza.getProperty("twoBitPath")));
        config.setNameSet("ucsc");
        config.setWholeGenomeView(false);

        if (this.genomeStanza.hasProperty("defaultPos")) {
            config.setDefaultPos(this.genomeStanza.getProperty("defaultPos"));
        }

        config.setDescription(config.getId());

        if (this.genomeStanza.hasProperty("blat")) {
            config.setBlat(getDataURL(this.genomeStanza.getProperty("blat")));
        }
        if (this.genomeStanza.hasProperty("chromAliasBb")) {
            config.setChromAliasBbURL(getDataURL(this.genomeStanza.getProperty("chromAliasBb")));
        }
        if (this.genomeStanza.hasProperty("twoBitBptURL")) {
            config.setTwoBitBptURL(getDataURL(this.genomeStanza.getProperty("twoBitBptURL")));
        }

        if (this.genomeStanza.hasProperty("twoBitBptUrl")) {
            config.setTwoBitBptURL(getDataURL(this.genomeStanza.getProperty("twoBitBptUrl")));
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
            config.setInfoURL(getDataURL(this.genomeStanza.getProperty("htmlPath")));
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
                config.setCytobandBbURL(getDataURL(t.getProperty("bigDataUrl")));
                break;
            }
        }

        // Tracks.  To prevent loading tracks set `includeTrackGroups`to false or "none"
        if (includeTracks) {
            Function<Stanza, Boolean> filter = (t) -> !Hub.filterTracks.contains(t.name) &&
                    (!"hide".equals(t.getProperty("visibility")));
            config.setTracks(this.getTrackConfigs(filter));
        }

        // config.trackConfigurations = this.#getGroupedTrackConfigurations()

        return config;
    }

    private String getDataURL(String relativeURL) {
        return relativeURL.startsWith("/") ? this.host + relativeURL : this.baseURL + relativeURL;
    }

    public List<TrackConfigGroup> getGroupedTrackConfigurations() {

        // Build map of group objects
        Map<String, TrackConfigGroup> groupMap = new HashMap<>();
        groupMap.put("other", new TrackConfigGroup("other", "Other", Integer.MAX_VALUE, false));
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

        final List<TrackConfig> trackConfigs = this.getTrackConfigs(stanza -> !stanza.name.equals("cytoBandIdeo"));
        for (TrackConfig c : trackConfigs) {
            String groupName = c.getGroup() != null ? c.getGroup() : "other";
            final TrackConfigGroup trackConfigGroup = groupMap.get(groupName);
            int priority = trackConfigGroup.priority;
            trackConfigGroup.tracks.add(c);

            String parentName = c.getStanzaParent();
            if (parentName != null && trackStanzaMap.containsKey(parentName)) {
                // Create a contingent container, will be used if a sufficient # of tracks belong to this container
                TrackConfigGroup container = parentCache.get(parentName);
                if (container == null) {
                    Stanza s = trackStanzaMap.get(parentName);
                    String label = trackConfigGroup.label + " - " + s.getProperty("shortLabel");
                    container = new TrackConfigGroup(parentName, label, priority + 1, false);
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


    TrackConfig getTrackConfig(Stanza t) {

        String format = t.format();
        String url = getDataURL(t.getProperty("bigDataUrl"));
        TrackConfig config = new TrackConfig(url);

        config.setPanelName(IGV.DATA_PANEL_NAME);

        config.setId(t.getProperty("track"));
        config.setName(t.getProperty("shortLabel"));

        // TODO -- work on recognizing big* formats
        // config.format = t.format();

        config.setUrl(getDataURL(t.getProperty("bigDataUrl")));

        // Expanded display mode does not work well in IGV desktop for some tracks
        //config.displayMode = t.displayMode();

        if ("vcfTabix".equals(format)) {
            config.setIndexURL(config.getUrl() + ".tbi");
        }

        if (t.hasProperty("longLabel") && t.hasProperty("html")) {
            config.setDescription("<a target=\"_blank\" href=\"" + (getDataURL(t.getProperty("html")) + "\">" + t.getProperty("longLabel") + "</a>"));
        } else if (t.hasProperty("longLabel")) {
            config.setDescription(t.getProperty("longLabel"));
        }

        if (t.hasProperty("html")) {
            config.setHtml(getDataURL(t.getProperty("html")));
        }

        String vizProperty = t.getProperty("visibility");

        if (vizProperty != null && vizModeMap.containsKey(vizProperty)) {
            config.setDisplayMode(vizModeMap.get(vizProperty));
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
            config.setTrixURL(getDataURL(t.getProperty("searchTrix")));
        }

        if (t.hasProperty("group")) {
            config.setGroup(t.getProperty("group"));
        }

        if (t.parent != null) {
            config.setStanzaParent(t.parent.name);
        }

        return config;
    }

    public String getUrl() {
        return url;
    }


    static class Stanza {

        private static Set<String> parentOverrideProperties = new HashSet<>(Arrays.asList("visibility", "priority", "group"));
        private final String type;
        private final String name;
        private Stanza parent;
        private Map<String, String> properties;


        Stanza(String type, String name) {

            this.type = type;
            this.name = name;
            this.properties = new HashMap<>();
        }

        public String getType() {
            return type;
        }

        public String getName() {
            return name;
        }

        void setProperty(String key, String value) {
            this.properties.put(key, value);
        }

        String getProperty(String key) {

            if (parentOverrideProperties.contains(key) && this.parent != null && this.parent.hasProperty(key)) {
                return this.parent.getProperty(key);
            } else if (this.properties.containsKey(key)) {
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

        public Stanza getParent() {
            return parent;
        }

        public void setParent(Stanza parent) {
            this.parent = parent;
        }
    }

    static String firstWord(String str) {
        return Globals.whitespacePattern.split(str)[0];
    }

}
