/*
 * Created by JFormDesigner on Thu Oct 31 22:31:02 EDT 2013
 */

package org.igv.encode;

import org.igv.Globals;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.ui.IGV;
import org.igv.ui.action.BrowseEncodeAction;
import org.igv.util.Pair;
import org.igv.util.ParsingUtils;

import javax.swing.text.NumberFormatter;
import java.awt.*;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Jim Robinson
 */
public class EncodeTrackChooserFactory {

    private static Logger log = LogManager.getLogger(EncodeTrackChooserFactory.class);

    private static Map<String, TrackChooser> instanceMap = Collections.synchronizedMap(new HashMap<>());
    private static NumberFormatter numberFormatter = new NumberFormatter();

    private static String ENCODE_HOST = "https://www.encodeproject.org";
    private static Set<String> filteredColumns = new HashSet(Arrays.asList(
            "ID", "Assembly", "HREF", "path",
            "url", "Project", "name", "color", "altColor"));

    private static List<String> filteredExtensions = Arrays.asList("tsv", "tsv.gz");

    private static Map<String, String> speciesNames = Map.of(
            "ce10", "Caenorhabditis elegans",
            "ce11", "Caenorhabditis elegans",
            "dm3", "Drosophila melanogaster",
            "dm6", "Drosophila melanogaster",
            "GRCh38", "Homo sapiens",
            "hg19", "Homo sapiens",
            "mm10", "Mus musculus",
            "mm9", "Mus musculus"
    );

    static HashSet<String> ucscSupportedGenomes = new HashSet<>(Arrays.asList("hg19", "mm9"));
    static HashSet<String> supportedGenomes = new HashSet<>(
            Arrays.asList("ce10", "ce11", "dm3", "dm6", "GRCh38", "hg19", "mm10", "mm9"));
    static HashSet<String> hicGenomes = new HashSet<>(
            Arrays.asList("GRCh38", "hg19", "mm10", "mm9"));

    /**
     * Return a new or cached instance of a track chooser for the given genome and type.
     *
     * @param genomeId
     * @param type
     * @return
     * @throws IOException
     */
    public synchronized static TrackChooser getInstance(String genomeId, BrowseEncodeAction.Type type) throws IOException {

        String encodeGenomeId = getEncodeGenomeID(genomeId);
        String key = encodeGenomeId + type.toString();
        TrackChooser instance = instanceMap.get(key);
        if (instance == null) {
            Pair<List<String>, List<FileRecord>> records = getEncodeFileRecords(encodeGenomeId, type);
            if (records == null) {
                return null;
            }
            Frame parent = IGV.hasInstance() ? IGV.getInstance().getMainFrame() : null;
            final List<String> headings = records.getFirst();
            final List<FileRecord> rows = records.getSecond();
            final String title = getDialogTitle(genomeId, type);
            instance = new TrackChooser(parent, headings, rows, title);
            instanceMap.put(key, instance);
        }
        return instance;
    }

    private static String getDialogTitle(String genomeId, BrowseEncodeAction.Type type) {

        if (type == BrowseEncodeAction.Type.UCSC) {
            return "ENCODE data hosted at UCSC (2012)";
        } else if (type == BrowseEncodeAction.Type.FOUR_DN) {
            return "4DN";
        } else {
            switch (type) {
                case SIGNALS_CHIP:
                    return "ENCODE CHiP Seq - Signals";
                case SIGNALS_OTHER:
                    return "ENCODE Non CHiP Data - Signals";
                default:
                    return "ENCODE";
            }
        }
    }

    public static boolean genomeSupportedUCSC(String genomeId) {
        return genomeId != null && ucscSupportedGenomes.contains(getEncodeGenomeID(genomeId));
    }

    public static boolean hicSupportedUCSC(String genomeId) {
        return genomeId != null && hicGenomes.contains(getEncodeGenomeID(genomeId));
    }

    public static boolean genomeSupported(String genomeId) {
        return genomeId != null && supportedGenomes.contains(getEncodeGenomeID(genomeId));
    }


    private static String getEncodeGenomeID(String genomeId) {
        switch (genomeId) {
            case "hg38":
            case "hg38_1kg":
                return "GRCh38";
            case "b37":
            case "1kg_v37":
                return "hg19";
            default:
                return genomeId;
        }

    }

    private static Pair<List<String>, List<FileRecord>> getEncodeFileRecords(String genomeId, BrowseEncodeAction.Type type) throws IOException {

        try (InputStream is = getStreamFor(genomeId, type)) {
            if (is == null) {
                return null;
            }
            return parseRecords(is, type, genomeId);
        }
    }

    private static InputStream getStreamFor(String genomeId, BrowseEncodeAction.Type type) throws IOException {
        if (type == BrowseEncodeAction.Type.UCSC) {
            return EncodeTrackChooserFactory.class.getResourceAsStream("encode." + genomeId + ".txt");
        } else if (type == BrowseEncodeAction.Type.FOUR_DN) {
            String url = PreferencesManager.getPreferences().get(Constants.FOUR_DN_FILELIST_URL) + "4dn_" + genomeId + "_tracks.txt";
            return ParsingUtils.openInputStream(url);
        } else {
            String root = PreferencesManager.getPreferences().get(Constants.ENCODE_FILELIST_URL) + genomeId + ".";
            String url = null;
            switch (type) {
                case SIGNALS_CHIP:
                    url = root + "signals.chip.txt.gz";
                    break;
                case SIGNALS_OTHER:
                    url = root + "signals.other.txt.gz";
                    break;
                case HIC:
                    url = root + "hic.txt.gz";
                    break;
                case OTHER:
                    url = root + "other.txt.gz";
                    break;
            }
            if (url == null) {
                throw new RuntimeException("Unknown encode data collection type: " + type);
            }
            return ParsingUtils.openInputStream(url);
        }
    }

    private static Pair parseRecords(InputStream is, BrowseEncodeAction.Type type, String genomeId) throws IOException {

        BufferedReader reader = new BufferedReader(new InputStreamReader(is));

        String[] headers = Globals.tabPattern.split(reader.readLine());

        int pathColumn = switch (type) {
            case UCSC, FOUR_DN -> 0;
            default -> Arrays.asList(headers).indexOf("HREF");
        };

        List<FileRecord> records = new ArrayList<>(20000);
        String nextLine;
        while ((nextLine = reader.readLine()) != null) {
            if (!nextLine.startsWith("#")) {

                String[] tokens = Globals.tabPattern.split(nextLine, -1);
                String path = switch (type) {
                    case UCSC, FOUR_DN -> tokens[pathColumn];
                    default -> ENCODE_HOST + tokens[pathColumn];
                };

                if (filteredExtensions.stream().anyMatch(e -> path.endsWith(e))) {
                    continue;
                }

                Map<String, String> attributes = new HashMap<>();
                for (int i = 0; i < headers.length; i++) {
                    String value = i < tokens.length ? tokens[i] : "";
                    if (value.length() > 0) {
                        attributes.put(headers[i], shortenField(value, genomeId));
                    }
                }
                final FileRecord record = new FileRecord(path, attributes);
                records.add(record);
            }
        }

        List<String> filteredHeaders = Arrays.stream(headers).filter(h -> !filteredColumns.contains(h)).collect(Collectors.toList());

        return new Pair(filteredHeaders, records);
    }

    private static String shortenField(String value, String genomeId) {
        String species = speciesNames.get(genomeId);
        return species == null ?
                value :
                value.replace("(" + species + ")", "").replace(species, "").trim();
    }


    /**
     * Main function for testing only
     *
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {
        getInstance("hg19", BrowseEncodeAction.Type.UCSC).setVisible(true);
    }

}
