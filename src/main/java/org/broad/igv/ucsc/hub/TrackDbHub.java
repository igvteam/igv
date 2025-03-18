package org.broad.igv.ucsc.hub;

import org.broad.igv.feature.genome.load.TrackConfig;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.StringUtils;

import java.util.*;
import java.util.stream.Collectors;

import static org.broad.igv.ucsc.hub.Hub.getPriority;

public class TrackDbHub {


    static Set supportedTypes = new HashSet(Arrays.asList("bigbed", "bigwig", "biggenepred", "vcftabix", "refgene",
            "bam", "sampleinfo", "vcf.list", "ucscsnp", "bed", "tdf", "gff", "gff3", "gtf", "vcf", "vcfphasedtrio",
            "bigdbsnp"));

    static Set filterTracks = new HashSet(Arrays.asList("cytoBandIdeo", "assembly", "gap", "gapOverlap", "allGaps",
            "cpgIslandExtUnmasked", "windowMasker"));

    static Map<String, String> vizModeMap = Map.of(
            "pack", "EXPANDED",
            "full", "EXPANDED",
            "squish", "SQUISHED",
            "dense", "COLLAPSED");

    static Map<String, String> typeFormatMap = Map.of(
            "vcftabix", "vcf",
            "vcfphasedtrio", "vcf",
            "bigdbsnp", "bigbed"
    );


    List<Stanza> trackStanzas;
    List<Stanza> groupStanzas;
    List<TrackConfigContainer> groupTrackConfigs;

    public TrackDbHub(List<Stanza> trackStanzas, List<Stanza> groupStanzas) {
        this.groupStanzas = groupStanzas;
        this.trackStanzas = trackStanzas;
    }

    public List<TrackConfigContainer> getGroupedTrackConfigurations(String hubName) {

        if (groupTrackConfigs == null) {
            // Build map of group objects
            groupTrackConfigs = new ArrayList<>();
            Map<String, TrackConfigContainer> trackContainers = new HashMap<>();

            // Create a group for tracks without a group
            TrackConfigContainer nullContainer = new TrackConfigContainer(hubName, hubName, 0, true);
            groupTrackConfigs.add(nullContainer);

            if (this.groupStanzas != null) {
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

                if("TFrPeakClusters".equals(s.getProperty("track"))) {
                    System.out.println("TFrPeakClusters");
                }


                boolean isContainer = (s.hasOwnProperty("superTrack") && !s.hasOwnProperty("bigDataUrl")) ||
                        s.hasOwnProperty("compositeTrack") ||
                        s.hasOwnProperty("view") ||
                        (s.hasOwnProperty("container") && s.getOwnProperty("container").equals("multiWig"));

                // Find parent, if any.  "group" containers can be implicit, all other types should be explicitly
                // defined before their children
                TrackConfigContainer parent = null;

                if (s.hasOwnProperty("parent")) {
                    parent =  trackContainers.get(s.getOwnProperty("parent"));
                }

                if (parent == null && s.hasProperty("group")) {
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
                    String label = s.hasProperty("longLabel") ? s.getProperty("longLabel") : s.getProperty("shortLabel");
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
                        supportedTypes.contains(s.type().toLowerCase())) {

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

        String url = t.getProperty("bigDataUrl");
        TrackConfig config = new TrackConfig(url);

        String type = t.type();
        if (type != null) {
            String tlc = type.toLowerCase();
            String format = typeFormatMap.containsKey(tlc) ? typeFormatMap.get(tlc) : type;
            config.setFormat(format.toLowerCase());
        }

        config.setId(t.getProperty("track"));
        config.setName(t.getProperty("shortLabel"));

        // A rather complex workaround for some composite tracks
        String longLabel = t.getOwnProperty("longLabel");
        if (longLabel == null) {
            String inheritedLongLabel = t.getProperty("longLabel");
            longLabel = inheritedLongLabel != null && inheritedLongLabel.length() > config.getName().length() ? inheritedLongLabel : config.getName();
        }
        config.setLongLabel(longLabel);

        // TODO -- work on recognizing big* formats
        // config.format = t.format();

        config.setUrl(t.getProperty("bigDataUrl"));

        if (t.hasProperty("bigDataIndex")) {
            config.setIndexURL(t.getProperty("bigDataIndex"));
        } else if (t.type().equals("vcfTabix")) {
            config.setIndexURL(t.getProperty("bigDataUrl") + ".tbi");
        }

        // Expanded display mode does not work well in IGV desktop for some tracks
        //config.displayMode = t.displayMode();

        if (t.hasProperty("longLabel")) {
            config.setDescription(t.getProperty("longLabel"));
        }

        if (t.hasProperty("html")) {
            config.setHtml(t.getProperty("html"));
        }

        String vizProperty = t.getProperty("visibility");
        if (vizProperty != null && vizModeMap.containsKey(vizProperty)) {
            config.setDisplayMode(vizModeMap.get(vizProperty));
        }
        config.setVisible(vizProperty != null && !("hide".equals(vizProperty)));

        if (t.hasProperty("maxWindowToDraw")) {
            long maxWindow = Long.parseLong(t.getProperty("maxWindowToDraw"));
            int vizWindow = Math.min(Integer.MAX_VALUE, (int) maxWindow);
            config.setVisibilityWindow(vizWindow);
        }

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
        if (t.hasProperty("metadata")) {
            Map<String, String> metadata = parseMetadata(t.getProperty("metadata"));
            config.setAttributes(metadata);
        }
        if(t.hasProperty("defaultLabelFields")) {
            config.setLabelField(t.getProperty("defaultLabelFields").split(",")[0]);
        } else if (t.hasProperty("labelFields")) {
            config.setLabelField(t.getProperty("labelFields").split(",")[0]);
        }

        if (t.parent != null) {
            config.setStanzaParent(t.getAncestor().name);
        }

        return config;
    }

    //metadata differentiation="10 hour" treatment=X donor=A lab="List Meta Lab" data_set_id=ucscTest1 access=group assay=long-RNA-seq enriched_in=exon life_stage=postpartum species="Homo sapiens" ucsc_db=hg38
    static Map<String, String> parseMetadata(String metadata) {

        Map<String, String> attrs = new HashMap();
        while (metadata.length() > 0) {

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
}
