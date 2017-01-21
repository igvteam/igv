package org.broad.igv.prefs;

import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.ui.AboutDialog;
import org.broad.igv.ui.IGVCommandBar;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.event.IGVEventBus;
import org.broad.igv.ui.event.IGVEventObserver;
import org.broad.igv.util.ParsingUtils;

import java.awt.*;
import java.io.*;
import java.util.*;

import static org.broad.igv.prefs.Constants.*;
import static org.broad.igv.prefs.Constants.SESSION_RELATIVE_PATH;
import static org.broad.igv.prefs.Constants.SHOW_LOS;

/**
 * Created by jrobinso on 1/21/17.
 */
public class PreferencesManager implements IGVEventObserver {

    private static Logger log = Logger.getLogger(PreferencesManager.class);

    private static IGVPreferences genericPreferences;

    private static Map<String, IGVPreferences> preferencesMap = Collections.synchronizedMap(new HashMap<>());

    private static String prefFile;  // User preferences file


    private PreferencesManager() {
    }

    private static PreferencesManager theInstance = new PreferencesManager();

    public static IGVPreferences getPreferences(String category) {

        if (genericPreferences == null) {
            loadUserPreferences();
            IGVEventBus.getInstance().subscribe(PreferencesChangeEvent.class, theInstance);
        }
        if (preferencesMap.containsKey(category)) {
            return preferencesMap.get(category);
        } else {
            return genericPreferences;
        }
    }

    public static IGVPreferences getPreferences() {
        if (genericPreferences == null) {
            loadUserPreferences();
            IGVEventBus.getInstance().subscribe(PreferencesChangeEvent.class, theInstance);
        }
        return genericPreferences;
    }


    public static void setPrefsFile(String prefsFile) {
        prefFile = prefsFile;
    }


    private static void loadUserPreferences() {

        try {
            if (prefFile == null) {
                prefFile = DirectoryManager.getPreferencesFile().getAbsolutePath();
            }

            Map<String, Map<String, String>> userPrefs = load(prefFile);

            IGVPreferences defaultPreferences = new IGVPreferences(genericDefaults, null);
            genericPreferences = new IGVPreferences(userPrefs.get(NULL_CATEGORY), defaultPreferences);

            // If no user preferences are set for 3rd gen and RNA seed with reasonable defaults
            Map<String, String> thirdGenPreferences = userPrefs.get(THIRD_GEN);
            if(thirdGenPreferences == null || thirdGenPreferences.isEmpty()) {
                thirdGenPreferences = thirdGenDefaults;
            }
            preferencesMap.put(THIRD_GEN, new IGVPreferences(thirdGenPreferences, genericPreferences));

            Map<String, String> rnaPreferences = userPrefs.get(RNA);
            if(rnaPreferences == null || rnaPreferences.isEmpty()) {
                rnaPreferences = rnaDefaults;
            }
            preferencesMap.put(RNA, new IGVPreferences(rnaPreferences, genericPreferences));

        } catch (Exception e) {
            log.error("Error loading preferences file: " + prefFile, e);
            genericPreferences = new IGVPreferences();
        }

    }


    public static void loadOverrides(String overridePropertyFilePath) {

        if (genericPreferences == null) {
            loadUserPreferences();
        }

        Map<String, Map<String, String>> overrides = load(overridePropertyFilePath);

        for(Map.Entry<String, Map<String, String>> entry : overrides.entrySet()) {

            IGVPreferences prefs = preferencesMap.containsKey(entry.getKey()) ?
                    preferencesMap.get(entry.getKey()) :
                    genericPreferences;

            prefs.addOverrides(entry.getValue());
        }

    }

    private static Map<String, Map<String, String>> load(String prefFileName) {

        Map<String, Map<String, String>> prefMap = new HashMap<>();
        BufferedReader reader = null;
        try {
            reader = ParsingUtils.openBufferedReader(prefFileName);
            String nextLine = null;
            String category = NULL_CATEGORY;
            while ((nextLine = reader.readLine()) != null) {
                if (nextLine.startsWith("##")) {
                    category = nextLine.substring(2).trim();
                } else {
                    Map<String, String> prefs = prefMap.get(category);
                    if (prefs == null) {
                        prefs = Collections.synchronizedMap(new HashMap<>());
                        prefMap.put(category, prefs);
                    }

                    int idx = nextLine.indexOf('=');
                    if (idx > 0) {
                        String key = nextLine.substring(0, idx);
                        if (aliasTable.containsKey(key)) {
                            key = aliasTable.get(key);
                        }
                        String value = nextLine.substring(idx + 1);
                        prefs.put(key, value);
                    }
                }
            }
        } catch (IOException e) {
            log.error("Error loading preferences", e);
        } finally {
            try {
                if (reader != null) {
                    reader.close();
                }
            } catch (IOException ex) {
                log.error("Error closing preferences file", ex);
            }
        }
        return prefMap;
    }

    private synchronized void storePreferences() {

        if (genericPreferences != null) {
            FileWriter fileWriter = null;
            try {
                fileWriter = new FileWriter(prefFile);
                PrintWriter pw = new PrintWriter(new BufferedWriter(fileWriter));

                genericPreferences.print(pw);

                for(Map.Entry<String, IGVPreferences> entry : preferencesMap.entrySet()) {

                    pw.println();
                    pw.println("##" + entry.getKey());
                    entry.getValue().print(pw);

                }

                pw.flush();
                pw.close();
            } catch (IOException e) {
                log.error("Error loading preferences", e);
            } finally {

                if (fileWriter != null) {
                    try {
                        fileWriter.close();
                    } catch (IOException e) {
                        // Ignore
                    }
                }
            }

        }
    }


    static Hashtable<String, String> aliasTable = new Hashtable<String, String>();

    static {
        aliasTable.put("SAM>SORT_OPTION", "SAM.SORT_OPTION");
        aliasTable.put("FLAKING_REGIONS", "FLANKING_REGION");
    }

    @Override
    public void receiveEvent(Object event) {
        if (event instanceof PreferencesChangeEvent) {
            storePreferences();
        }
    }


    static Map<String, String> genericDefaults = new HashMap<>();
    static Map<String, String> thirdGenDefaults = new HashMap<>();
    static Map<String, String> rnaDefaults = new HashMap<>();

    static {

        genericDefaults.put(MUTATION_INDEL_COLOR_KEY, "0,200,0");
        genericDefaults.put(MUTATION_MISSENSE_COLOR_KEY, "170,20,240");
        genericDefaults.put(MUTATION_NONSENSE_COLOR_KEY, "50,30,75");
        genericDefaults.put(MUTATION_SPLICE_SITE_COLOR_KEY, "150,0,150");
        genericDefaults.put(MUTATION_SYNONYMOUS_COLOR_KEY, "200,170,200");
        genericDefaults.put(MUTATION_TARGETED_REGION_COLOR_KEY, "236,155,43");
        genericDefaults.put(MUTATION_UNKNOWN_COLOR_KEY, "0,180,225");
        //     * Nico's labels:   Truncating, Non-coding_Transcript, Other_AA_changing, Other_likely_neutral.
        genericDefaults.put("MUTATION_Truncating_COLOR", "150,0,0");
        genericDefaults.put("MUTATION_Non-coding_Transcript_COLOR", "0,0,150");
        genericDefaults.put("MUTATION_Other_AA_changing_COLOR", "0,150,150");
        genericDefaults.put("MUTATION_Other_likely_neutral_COLOR", "225,180,225");


        genericDefaults.put(PROBE_MAPPING_KEY, "false");
        genericDefaults.put(PROBE_MAPPING_FILE, null);
        genericDefaults.put(USE_PROBE_MAPPING_FILE, "false");

        genericDefaults.put(SHOW_REGION_BARS, "false");

        genericDefaults.put(OVERLAY_MUTATION_TRACKS, "true");
        genericDefaults.put(SHOW_ORPHANED_MUTATIONS, "true");
        genericDefaults.put(COLOR_MUTATIONS, "false");
        genericDefaults.put(OVERLAY_MUTATIONS_WHOLE_GENOME, "true");
        genericDefaults.put(SHOW_SINGLE_TRACK_PANE_KEY, "false");
        genericDefaults.put(PORT_ENABLED, "true");
        genericDefaults.put(EXPAND_FEAUTRE_TRACKS, "false");
        genericDefaults.put(SHOW_ATTRIBUTE_VIEWS_KEY, "true");
        genericDefaults.put(SHOW_SINGLE_TRACK_PANE_KEY, "false");
        genericDefaults.put(SHOW_DEFAULT_TRACK_ATTRIBUTES, "false");

        genericDefaults.put(CHART_DRAW_TOP_BORDER, "false");
        genericDefaults.put(CHART_DRAW_BOTTOM_BORDER, "false");
        genericDefaults.put(CHART_COLOR_BORDERS, "true");
        genericDefaults.put(CHART_DRAW_TRACK_NAME, "false");
        genericDefaults.put(CHART_DRAW_Y_AXIS, "false");
        genericDefaults.put(CHART_AUTOSCALE, "false");
        genericDefaults.put(CHART_SHOW_DATA_RANGE, "true");
        genericDefaults.put(CHART_COLOR_TRACK_NAME, "true");
        genericDefaults.put(CHART_TRACK_HEIGHT_KEY, "40");
        genericDefaults.put(CHART_SHOW_ALL_HEATMAP, "false");

        genericDefaults.put(SAM_SHOW_DUPLICATES, "false");
        genericDefaults.put(SAM_QUICK_CONSENSUS_MODE, "false");
        genericDefaults.put(SAM_SHOW_SOFT_CLIPPED, "false");
        genericDefaults.put(SAM_FLAG_UNMAPPED_PAIR, "false");
        genericDefaults.put(SAM_AUTO_SORT, "false");
        genericDefaults.put(SAM_SHADE_CENTER, "true");
        genericDefaults.put(SAM_SHOW_REF_SEQ, "false");
        genericDefaults.put(SAM_SHOW_CENTER_LINE, "true");
        genericDefaults.put(SAM_SHOW_COV_TRACK, "true");
        genericDefaults.put(SAM_SHADE_BASES, AlignmentTrack.ShadeBasesOption.QUALITY.toString());
        genericDefaults.put(SAM_FILTER_ALIGNMENTS, "false");
        genericDefaults.put(SAM_FILTER_SECONDARY_ALIGNMENTS, "false");
        genericDefaults.put(SAM_FILTER_SUPPLEMENTARY_ALIGNMENTS, "false");
        genericDefaults.put(SAM_FILTER_FAILED_READS, "true");
        genericDefaults.put(SAM_DOWNSAMPLE_READS, "true");
        genericDefaults.put(SAM_SAMPLING_WINDOW, "50");
        genericDefaults.put(SAM_SAMPLING_COUNT, "100");
        genericDefaults.put(SAM_BASE_QUALITY_MIN, "5");
        genericDefaults.put(SAM_BASE_QUALITY_MAX, "20");
        genericDefaults.put(SAM_FILTER_URL, null);
        genericDefaults.put(SAM_HIDDEN_TAGS, "SA,MD,XA,RG");
        genericDefaults.put(SAM_QUALITY_THRESHOLD, "0");
        genericDefaults.put(SAM_ALLELE_THRESHOLD, "0.2f");
        genericDefaults.put(SAM_ALLELE_USE_QUALITY, "true");
        genericDefaults.put(SAM_MIN_INSERT_SIZE_THRESHOLD, "50");
        genericDefaults.put(SAM_MAX_INSERT_SIZE_THRESHOLD, "1000");
        genericDefaults.put(SAM_MIN_INSERT_SIZE_PERCENTILE, "0.5");
        genericDefaults.put(SAM_MAX_INSERT_SIZE_PERCENTILE, "99.5");
        genericDefaults.put(SAM_MAX_VISIBLE_RANGE, "30");
        genericDefaults.put(SAM_COLOR_BY, "UNEXPECTED_PAIR");
        genericDefaults.put(SAM_COLOR_BY_TAG, "");
        genericDefaults.put(SAM_GROUP_BY_TAG, "");
        genericDefaults.put(SAM_GROUP_BY_POS, "");
        genericDefaults.put(SAM_SORT_BY_TAG, "");
        genericDefaults.put(SAM_BISULFITE_CONTEXT, "CG");
        genericDefaults.put(SAM_COMPUTE_ISIZES, "true");
        genericDefaults.put(SAM_FLAG_ZERO_QUALITY, "true");
        genericDefaults.put(SAM_SHOW_JUNCTION_TRACK, "false");
        genericDefaults.put(SAM_JUNCTION_MIN_FLANKING_WIDTH, "0");
        genericDefaults.put(SAM_JUNCTION_MIN_COVERAGE, "1");
        genericDefaults.put(SAM_SHOW_JUNCTION_FLANKINGREGIONS, "true");
        genericDefaults.put(SAM_NOMESEQ_ENABLED, "false");
        genericDefaults.put(SAM_COUNT_DELETED_BASES_COVERED, "false");
        genericDefaults.put(SAM_FLAG_LARGE_INDELS, "true");
        genericDefaults.put(SAM_LARGE_INDELS_THRESHOLD, "1");
        genericDefaults.put(SAM_FLAG_CLIPPING, "false");
        genericDefaults.put(SAM_CLIPPING_THRESHOLD, "0");
        genericDefaults.put(SAM_SORT_OPTION, "NUCLEOTIDE");
        genericDefaults.put(SAM_GROUP_OPTION, "NONE");
        genericDefaults.put(SAM_SHOW_GROUP_SEPARATOR, "true");
        genericDefaults.put(SAM_COMPLETE_READS_ONLY, "false");
        genericDefaults.put(SAM_SHOW_ALL_BASES, "false");

        genericDefaults.put(SAM_REDUCED_MEMORY_MODE, "false");

        genericDefaults.put(SAM_HIDE_SMALL_INDEL_BP, "false");
        genericDefaults.put(SAM_SMALL_INDEL_BP_THRESHOLD, "0");

        genericDefaults.put(SAM_LINK_READS, "false");
        genericDefaults.put(SAM_LINK_TAG, "READNAME");

        genericDefaults.put(SAM_SHOW_ALIGNMENT_TRACK, "true");

        genericDefaults.put(BYPASS_FILE_AUTO_DISCOVERY, "false");

        genericDefaults.put(NORMALIZE_COVERAGE, "false");

        genericDefaults.put(SHOW_GENOME_SERVER_WARNING, "true");

        genericDefaults.put(SEARCH_ZOOM, "true");


        genericDefaults.put(GENOMES_SERVER_URL, Globals.DEFAULT_GENOME_URL);
        genericDefaults.put(OVERLAY_ATTRIBUTE_KEY, "LINKING_ID");
        genericDefaults.put(DEFAULT_GENOME, Globals.DEFAULT_GENOME);

        genericDefaults.put(USE_PROXY, "false");
        genericDefaults.put(PROXY_AUTHENTICATE, "false");
        genericDefaults.put(PORT_NUMBER, "60151");
        genericDefaults.put(TRACK_HEIGHT_KEY, "15");
        genericDefaults.put(FLANKING_REGION, "2000");

        genericDefaults.put(SHOW_SEQUENCE_TRANSLATION, "false");
        genericDefaults.put(MAX_SEQUENCE_RESOLUTION, "2");

        genericDefaults.put(AUTO_UPDATE_GENOMES, "true");

        genericDefaults.put(GWAS_TRACK_HEIGHT, "200");
        genericDefaults.put(GWAS_DESCRIPTION_CACHE_SIZE, "10000");
        genericDefaults.put(GWAS_MIN_POINT_SIZE, "3");
        genericDefaults.put(GWAS_MAX_POINT_SIZE, "7");
        genericDefaults.put(GWAS_USE_CHR_COLORS, "true");
        genericDefaults.put(GWAS_SINGLE_COLOR, "false");
        genericDefaults.put(GWAS_ALTERNATING_COLORS, "false");
        genericDefaults.put(GWAS_PRIMARY_COLOR, "69,101,183");
        genericDefaults.put(GWAS_SECONDARY_COLOR, "250,169,10");
        genericDefaults.put(GWAS_SHOW_AXIS, "true");

        genericDefaults.put(DEFAULT_FONT_SIZE, "10");
        genericDefaults.put(DEFAULT_FONT_FAMILY, "Arial");
        genericDefaults.put(DEFAULT_FONT_ATTRIBUTE, String.valueOf(Font.PLAIN));
        genericDefaults.put(SCALE_FONTS, "false");

        boolean isMac = System.getProperty("os.name").toLowerCase().startsWith("mac");
        genericDefaults.put(ENABLE_ANTIALISING, String.valueOf(isMac));

        genericDefaults.put(NAME_PANEL_WIDTH, "160");
        genericDefaults.put(BACKGROUND_COLOR, "250,250,250");

        genericDefaults.put(GENOME_SPACE_ENABLE, "true");
        genericDefaults.put(GENOME_SPACE_DM_SERVER, "https://dm.genomespace.org/datamanager/v1.0/");
        genericDefaults.put(GENOME_SPACE_ATM_SERVER, "https://atm.genomespace.org/atm/v1.0/");
        genericDefaults.put(GENOME_SPACE_IDENTITY_SERVER, "https://identitydev.genomespace.org:8444/identityServer/basic");

        genericDefaults.put(DB_ENABLED, "false");
        genericDefaults.put(DB_HOST, "");
        genericDefaults.put(DB_NAME, "");
        genericDefaults.put(DB_PORT, "-1");


        String defaultDataURL = Globals.DEFAULT_DATA_URL;
        Properties properties = new Properties();
        try {
            properties.load(AboutDialog.class.getResourceAsStream("/resources/about.properties"));
            String tmp = properties.getProperty("master-resource-url");
            if (tmp != null && !tmp.startsWith("@")) {
                defaultDataURL = tmp;
            }
        } catch (IOException e) {
            log.error("Error reading dataURL property", e);
        }

        genericDefaults.put(DATA_SERVER_URL_KEY, defaultDataURL);

        genericDefaults.put(CBIO_MUTATION_THRESHOLD, "1");
        genericDefaults.put(CBIO_AMPLIFICATION_THRESHOLD, "0.9");
        genericDefaults.put(CBIO_DELETION_THRESHOLD, "0.9");
        genericDefaults.put(CBIO_EXPRESSION_UP_THRESHOLD, "1.0");
        genericDefaults.put(CBIO_EXPRESSION_DOWN_THRESHOLD, "1.0");

        genericDefaults.put(TOOLTIP_INITIAL_DELAY, "50");
        genericDefaults.put(TOOLTIP_RESHOW_DELAY, "50");
        genericDefaults.put(TOOLTIP_DISMISS_DELAY, "60000");
        genericDefaults.put(DETAILS_BEHAVIOR_KEY, IGVCommandBar.SHOW_DETAILS_BEHAVIOR.HOVER.name());

        genericDefaults.put(SHOW_SIZE_WARNING, "true");

        genericDefaults.put(SKIP_VERSION, "");

        genericDefaults.put(COLOR_A, "0,150,0");
        genericDefaults.put(COLOR_C, "0,0,255");
        genericDefaults.put(COLOR_T, "255,0,0");
        genericDefaults.put(COLOR_G, "209,113,5");
        genericDefaults.put(COLOR_N, ColorUtilities.colorToString(Color.gray));
        genericDefaults.put(SAM_COLOR_A, "0,255,0");
        genericDefaults.put(SAM_COLOR_C, "0,0,255");
        genericDefaults.put(SAM_COLOR_T, "255,0,0");
        genericDefaults.put(SAM_COLOR_G, "209,113,5");
        genericDefaults.put(SAM_COLOR_N, ColorUtilities.colorToString(Color.gray.brighter()));

        genericDefaults.put(HOMREF_COLOR, "235,235,235");
        genericDefaults.put(HETVAR_COLOR, "0,0,255");
        genericDefaults.put(HOMVAR_COLOR, "0,245,255");
        genericDefaults.put(NOCALL_COLOR, "255,255,255");
        genericDefaults.put(AF_REF_COLOR, "0,0,220");
        genericDefaults.put(AF_VAR_COLOR, "255,0,0");

        genericDefaults.put(VARIANT_COLOR_BY_ALLELE_FREQ, "true");

        genericDefaults.put(SASHIMI_SHOW_COVERAGE, "true");

        genericDefaults.put(ENABLE_GOOGLE_MENU, "false");
        genericDefaults.put(SAVE_GOOGLE_CREDENTIALS, "true");

        genericDefaults.put(DEFAULT_VISIBILITY_WINDOW, "-1");

        genericDefaults.put(BLAT_URL, "http://genome.cse.ucsc.edu/cgi-bin/hgBlat");

        genericDefaults.put(GENE_LIST_BED_FORMAT, "false");

        genericDefaults.put(SESSION_RELATIVE_PATH, "false");

        genericDefaults.put(SHOW_LOS, "true");


        thirdGenDefaults.put(SAM_QUICK_CONSENSUS_MODE, "true");
        thirdGenDefaults.put(SAM_DOWNSAMPLE_READS, "false");
        thirdGenDefaults.put(SAM_MAX_VISIBLE_RANGE, "1000");
        thirdGenDefaults.put(SAM_FLAG_LARGE_INDELS, "true");
        thirdGenDefaults.put(SAM_LARGE_INDELS_THRESHOLD, "1");
        thirdGenDefaults.put(SAM_FLAG_CLIPPING, "true");
        thirdGenDefaults.put(SAM_CLIPPING_THRESHOLD, "0");


        rnaDefaults.put(SAM_MAX_VISIBLE_RANGE, "300");
        rnaDefaults.put(SAM_SHOW_JUNCTION_TRACK, "true");


    }
}
