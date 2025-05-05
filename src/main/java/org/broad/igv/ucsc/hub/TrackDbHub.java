package org.broad.igv.ucsc.hub;

import org.broad.igv.feature.genome.load.TrackConfig;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.util.StringUtils;

import java.util.*;
import java.util.stream.Collectors;

public class TrackDbHub {

    private static Logger log = LogManager.getLogger(TrackDbHub.class);

    static Set supportedTypes = new HashSet(Arrays.asList("bigbed", "bigwig", "biggenepred", "vcftabix", "refgene",
            "bam", "sampleinfo", "vcf.list", "ucscsnp", "bed", "tdf", "gff", "gff3", "gtf", "vcf", "vcfphasedtrio",
            "bigdbsnp", "rmask", "genepred", "wig", "bedgraph", "interact", "broadpeak", "narrowpeak", "gappedpeak",
            "gistic", "seg", "mut"));

    static Set filterTracks = new HashSet(Arrays.asList("cytoBandIdeo", "assembly"));

    static Map<String, String> vizModeMap = Map.of(
            "pack", "EXPANDED",
            "full", "EXPANDED",
            "squish", "SQUISHED",
            "dense", "COLLAPSED");

    static Map<String, String> typeFormatMap = Map.of(
            "vcftabix", "vcf",
            "vcfphasedtrio", "vcf",
            "bigdbsnp", "bigbed",
            "genepred", "refgene"
    );


    List<Stanza> trackStanzas;
    List<Stanza> groupStanzas;
    List<TrackConfigContainer> groupTrackConfigs;

    public TrackDbHub(List<Stanza> trackStanzas, List<Stanza> groupStanzas) {
        this.groupStanzas = groupStanzas;
        this.trackStanzas = trackStanzas;
    }

    public String findCytobandURL() {
        for (Stanza t : this.trackStanzas) {
            if ("cytoBandIdeo".equals(t.name) && t.hasProperty("bigDataUrl")) {
                return t.getProperty("bigDataUrl");
            }
        }
        return null;
    }

    public int getSupportedTrackCount() {
        int count = 0;
        for (Stanza t : this.trackStanzas) {
            if (!filterTracks.contains(t.name) &&
                    t.hasProperty("bigDataUrl") &&
                    t.format() != null &&
                    supportedTypes.contains(t.format().toLowerCase())) {
                count++;
            }
        }
        return count;
    }

    public List<TrackConfigContainer> getGroupedTrackConfigurations(String hubName) {

        if (groupTrackConfigs == null) {
            // Build map of group objects
            groupTrackConfigs = new ArrayList<>();
            Map<String, TrackConfigContainer> trackContainers = new HashMap<>();

            // Create a group for tracks without a group
            TrackConfigContainer nullContainer = new TrackConfigContainer(hubName, hubName, 0, true);
            groupTrackConfigs.add(nullContainer);

            boolean hasGroups = this.groupStanzas != null;
            if (hasGroups) {
                for (Stanza groupStanza : this.groupStanzas) {
                    String name = groupStanza.getProperty("name");
                    boolean defaultOpen = "0".equals(groupStanza.getProperty("defaultIsClosed"));
                    int priority = groupStanza.hasProperty("priority") ? getPriority(groupStanza.getProperty("priority")) : Integer.MAX_VALUE - 1;
                    final TrackConfigContainer container = new TrackConfigContainer(name, groupStanza.getProperty("label"), priority, defaultOpen);
                    trackContainers.put(name, container);
                    groupTrackConfigs.add(container);
                }
            }

            for (Stanza s : trackStanzas) {

                boolean isContainer = (s.hasOwnProperty("superTrack") && !s.hasOwnProperty("bigDataUrl")) ||
                        s.hasOwnProperty("compositeTrack") || s.hasOwnProperty("view") ||
                        (s.hasOwnProperty("container") && s.getOwnProperty("container").equals("multiWig"));

                // Find parent, if any.  "group" containers can be implicit, all other types should be explicitly
                // defined before their children
                TrackConfigContainer parent = null;

                if (s.hasOwnProperty("parent")) {
                    parent = trackContainers.get(s.getOwnProperty("parent"));
                }

                if (parent == null && hasGroups && s.hasProperty("group")) {
                    String groupName = s.getProperty("group");
                    if (trackContainers.containsKey(groupName)) {
                        parent = trackContainers.get(groupName);
                    } else {
                        TrackConfigContainer container = new TrackConfigContainer(groupName, groupName, 1000, true);
                        trackContainers.put(groupName, container);
                        groupTrackConfigs.add(container);
                        parent = container;
                    }
                }

                if (isContainer) {

                    String name = s.getProperty("track");
                    int priority = s.hasProperty("priority") ? getPriority(s.getProperty("priority")) : Integer.MAX_VALUE - 1;
                    boolean defaultOpen = "0".equals(s.getProperty("defaultIsClosed"));
                    String longLabel = s.getProperty("longLabel");
                    String label = longLabel != null && longLabel.length() < 50 ? longLabel : s.getProperty("shortLabel");
                    final TrackConfigContainer container = new TrackConfigContainer(name, label, priority, defaultOpen);
                    if (trackContainers.containsKey(name)) {
                        throw new RuntimeException("Duplicate track container: " + name);
                    }
                    trackContainers.put(name, container);

                    if (parent == null) {
                        // No parent or a superTrack => promote to top level
                        groupTrackConfigs.add(container);
                    } else {
                        parent.children.add(container);
                    }

                } else if (!filterTracks.contains(s.name) &&
                        s.hasProperty("bigDataUrl") &&
                        s.format() != null &&
                        supportedTypes.contains(s.format().toLowerCase())) {

                    final TrackConfig trackConfig = getTrackConfig(s);
                    if (parent != null) {
                        parent.tracks.add(trackConfig);
                    } else {
                        nullContainer.tracks.add(trackConfig);
                    }
                }
            }

            // Filter empty groups and sort
            for (TrackConfigContainer c : groupTrackConfigs) {
                c.trim();
            }
            groupTrackConfigs = groupTrackConfigs.stream().filter(t -> !t.isEmpty()).collect(Collectors.toList());


            Collections.sort(groupTrackConfigs, Comparator.comparingInt(o -> o.priority));
        }
        return groupTrackConfigs;
    }

    private TrackConfig getTrackConfig(Stanza t) {
        TrackConfig config = new TrackConfig(t.getProperty("bigDataUrl"));

        String type = t.format();
        if (type != null) {
            String format = typeFormatMap.getOrDefault(type.toLowerCase(), type);
            config.format = format.toLowerCase();
        }

        config.id = t.getProperty("track");
        config.name = t.getProperty("shortLabel");

        String longLabel = t.getOwnProperty("longLabel");
        if (longLabel == null) {
            longLabel = Optional.ofNullable(t.getProperty("longLabel"))
                    .filter(label -> label.length() > config.name.length())
                    .orElse(config.name);
        }
        config.longLabel = longLabel;

        config.url = t.getProperty("bigDataUrl");
        config.indexURL = t.hasProperty("bigDataIndex") ? t.getProperty("bigDataIndex") :
                "vcfTabix".equals(t.format()) ? t.getProperty("bigDataUrl") + ".tbi" : null;

        config.description = t.getProperty("longLabel");
        config.html = t.getProperty("html");

        String vizProperty = t.getProperty("visibility");
        if (vizProperty != null) {
            config.displayMode = vizModeMap.getOrDefault(vizProperty, null);
            config.visible = !"hide".equals(vizProperty);
        }

        if (t.hasProperty("maxWindowToDraw")) {
            long maxWindowToDraw = Long.parseLong(t.getProperty("maxWindowToDraw"));
            if (maxWindowToDraw > Integer.MAX_VALUE) {
                maxWindowToDraw = Integer.MAX_VALUE;
            }
            config.visibilityWindow = (int) maxWindowToDraw;
        }

        if (t.hasProperty("autoScale")) {
            String value = t.getProperty("autoScale").toLowerCase();
            config.autoscale = "on".equals(value);
            if ("group".equals(value) && t.getParent() != null) {
                config.autoscaleGroup = t.getParent().name;
            }
        }

        if (t.hasProperty("maxHeightPixels")) {
            String[] tokens = t.getProperty("maxHeightPixels").split(":");
            config.maxHeight = Integer.parseInt(tokens[0]);
            if (tokens.length > 1) config.height = Integer.parseInt(tokens[1]);
            if (tokens.length > 2) config.minHeight = Integer.parseInt(tokens[2]);
        }

        config.color = parseColor(t.getProperty("color"));
        config.altColor = parseColor(t.getProperty("altColor"));
        if (t.hasProperty("viewLimits")) {
            String[] tokens = t.getProperty("viewLimits").split(":");
            if (tokens.length == 1) {
                config.max = Float.parseFloat(tokens[0]);
            } else if (tokens.length > 1) {
                config.min = Float.parseFloat(tokens[0]);
                config.max = Float.parseFloat(tokens[1]);
            }
        }

        config.infoURL = t.getProperty("url");
        config.searchIndex = t.getProperty("searchIndex");
        config.trixURL = t.getProperty("searchTrix");
        config.group = t.getProperty("group");
        if(t.hasProperty("metadata")) {
            config.attributes = parseMetadata(t.getProperty("metadata"));
        }
        String labelFields = t.hasProperty("defaultLabelFields") ?
                t.getProperty("defaultLabelFields") :
                t.getProperty("labelFields");
        if (labelFields != null) {
            config.labelField = labelFields.split(",")[0];
        }

        if (t.getParent() != null) {
            config.stanzaParent = t.getAncestor().name;
        }

        return config;
    }

    private String parseColor(String color) {
        return color != null && color.contains(",") ? "rgb(" + color + ")" : color;
    }

    //metadata differentiation="10 hour" treatment=X donor=A lab="List Meta Lab" data_set_id=ucscTest1 access=group assay=long-RNA-seq enriched_in=exon life_stage=postpartum species="Homo sapiens" ucsc_db=hg38
    static Map<String, String> parseMetadata(String metadata) {

        Map<String, String> attrs = new HashMap();
        int lastMetadataLength = 0;
        while (metadata != null && metadata.length() > 0) {
            if (lastMetadataLength == metadata.length()) {
                break;  // prevent infinite loop
            }
            lastMetadataLength = metadata.length();
            try {
                int idx = metadata.indexOf("=");
                if (idx == -1 || idx == metadata.length() - 1) {
                    break;
                }
                int idx2;
                String key = StringUtils.stripQuotes(capitalize(metadata.substring(0, idx)));
                String value;

                if ('"' == metadata.charAt(idx + 1)) {
                    idx++;
                    idx2 = metadata.indexOf("\" ", idx + 1);
                    value = idx2 > 0 ? metadata.substring(idx + 1, idx2) : metadata.substring(idx + 1);
                    idx2++;
                } else {
                    idx2 = metadata.indexOf(" ", idx + 1);
                    if (idx2 == -1) {
                        idx2 = metadata.length();
                    }
                    value = metadata.substring(idx + 1, idx2);
                }
                value = StringUtils.stripQuotes(value);
                if (value.endsWith("\"")) {
                    value = value.substring(0, value.length() - 1);
                }
                if (value.startsWith("<") && value.endsWith(">")) {
                    value = htmlText(value);
                }
                attrs.put(key, value);
                if (idx2 == metadata.length()) {
                    break;
                }
                metadata = idx2 > 0 ? metadata.substring(idx2 + 1).trim() : "";
            } catch (Exception e) {
                // We don't want to fail parsing the hub due to a failure parsing metadata.  Also we don't want to
                // overwhelm the log.  Metatdata is or marginal importance in IGV.
            }
        }
        return attrs;
    }

    private static String capitalize(final String line) {
        return Character.toUpperCase(line.charAt(0)) + line.substring(1);
    }


    static String htmlText(String html) {
        // Assumes a pattern like <span style="color:#C58DAA">Digestive</span>
        int idx1 = html.indexOf('>');
        int idx2 = html.indexOf('<', idx1);
        if (idx1 > 0 && idx2 > idx1) {
            return html.substring(idx1 + 1, idx2);
        } else {
            return html;
        }
    }

    /**
     * Return the priority for the group.  The priority format is uncertain, but extends to at least 2 levels (e.g. 3.4).
     * Ignore levels > 3
     *
     * @param priorityString Priority as a string (e.g. 3.4)
     * @return A priority as an integer
     */

    static int getPriority(String priorityString) {
        try {
            String[] tokens = priorityString.trim().split("\\.");
            int p = Integer.parseInt(tokens[0]) * 100;
            if (tokens.length > 1) {
                p += Integer.parseInt(tokens[1]) * 10;
            }
            if (tokens.length > 2) {
                p += Integer.parseInt(tokens[2]);
            }
            return p;
        } catch (Exception e) {
            log.error("Error parsing priority string: " + priorityString, e);
            return Integer.MAX_VALUE;
        }
    }

}
