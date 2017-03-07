package org.broad.igv.prefs;

import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.ui.AboutDialog;
import org.broad.igv.ui.IGVCommandBar;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.event.IGVEventObserver;
import org.broad.igv.util.ParsingUtils;

import java.awt.*;
import java.io.*;
import java.util.*;
import java.util.List;

import static org.broad.igv.prefs.Constants.*;
import static org.broad.igv.prefs.Constants.SESSION_RELATIVE_PATH;
import static org.broad.igv.prefs.Constants.SHOW_LOS;

/**
 * Created by jrobinso on 1/21/17.
 */
public class PreferencesManager implements IGVEventObserver {

    public static final String SEPARATOR_KEY = "---";
    public static final String INFO_KEY = "info";

    private static List<PreferenceGroup> preferenceGroupList;
    private static Logger log = Logger.getLogger(PreferencesManager.class);

    private static Map<String, IGVPreferences> preferencesMap = Collections.synchronizedMap(new HashMap<>());

    private static String prefFile;  // User preferences file

    static Hashtable<String, String> aliasTable = new Hashtable<String, String>();
    static {
        aliasTable.put("SAM>SORT_OPTION", "SAM.SORT_OPTION");
        aliasTable.put("FLAKING_REGIONS", "FLANKING_REGION");
    }


    private PreferencesManager() {

    }

    private static PreferencesManager theInstance = new PreferencesManager();

    public static IGVPreferences getPreferences(String category) {

        if (preferenceGroupList == null) {
            init();
        }
        if (preferencesMap.containsKey(category)) {
            return preferencesMap.get(category);
        } else {
            return preferencesMap.get(NULL_CATEGORY);
        }
    }

    private static void init() {

        try {
            preferenceGroupList = loadPreferenceList();

            Map<String, Map<String, String>> defaultPreferences = new HashMap<>();
            defaultPreferences.put(NULL_CATEGORY, new HashMap<>());
            defaultPreferences.put(RNA, new HashMap<>());
            defaultPreferences.put(THIRD_GEN, new HashMap<>());

            for (PreferenceGroup group : preferenceGroupList) {
                String category = group.category == null ? NULL_CATEGORY : group.category;
                Map<String, String> defs = defaultPreferences.get(category);
                if (defs == null) {
                    log.info("Unknown preference category: " + category);
                }
                for (Preference pref : group.preferences) {
                    defs.put(pref.getKey(), pref.getDefaultValue());
                }
            }

            IGVPreferences genericDefaults = new IGVPreferences(defaultPreferences.get(NULL_CATEGORY), null);
            IGVPreferences rnaDefaults = new IGVPreferences(defaultPreferences.get(RNA), genericDefaults);
            IGVPreferences thirdGenDefaults = new IGVPreferences(defaultPreferences.get(THIRD_GEN), genericDefaults);

            Map<String, Map<String, String>> userPrefs = loadUserPreferences();

            preferencesMap.put(NULL_CATEGORY, new IGVPreferences(userPrefs.get(NULL_CATEGORY), genericDefaults));
            preferencesMap.put(RNA, new IGVPreferences(userPrefs.get(RNA), rnaDefaults));
            preferencesMap.put(THIRD_GEN, new IGVPreferences(userPrefs.get(THIRD_GEN), thirdGenDefaults));


        } catch (IOException e) {
            e.printStackTrace();
        }

        loadUserPreferences();

        IGVEventBus.getInstance().subscribe(PreferencesChangeEvent.class, theInstance);
    }

    public static IGVPreferences getPreferences() {
        return getPreferences(NULL_CATEGORY);
    }

    public static Collection<IGVPreferences> getAllPreferences() {
        return preferencesMap.values();
    }

    public static void setPrefsFile(String prefsFile) {
        prefFile = prefsFile;
    }

    private static Map<String, Map<String, String>> loadUserPreferences() {

        try {
            if (prefFile == null) {
                prefFile = DirectoryManager.getPreferencesFile().getAbsolutePath();
            }
            return load(prefFile);
        } catch (Exception e) {
            log.error("Error loading preferences file: " + prefFile, e);
            return null;
        }

    }

    public static void updateAll(Map<String, Map<String, String>> preferenceMap) {
        for (Map.Entry<String, Map<String, String>> entry : preferenceMap.entrySet()) {
            IGVPreferences preferences = getPreferences(entry.getKey());
            if (preferences != null) {
                preferences.putAll(entry.getValue());
            }

        }
    }

    public static void loadOverrides(String overridePropertyFilePath) {

        if (preferencesMap.get(NULL_CATEGORY) == null) {
            loadUserPreferences();
        }

        Map<String, Map<String, String>> overrides = load(overridePropertyFilePath);

        for (Map.Entry<String, Map<String, String>> entry : overrides.entrySet()) {

            IGVPreferences prefs = preferencesMap.containsKey(entry.getKey()) ?
                    preferencesMap.get(entry.getKey()) :
                    preferencesMap.get(NULL_CATEGORY);

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

    public static List<PreferenceGroup> loadPreferenceList() throws IOException {

        List<PreferenceGroup> groupList = new ArrayList<>();
        List<Preference> prefList = null;

        BufferedReader reader = new BufferedReader(new InputStreamReader(PreferenceEditorFX.class.getResourceAsStream("/org/broad/igv/prefs/preferences.tab")));
        String nextLine;
        String group = null;

        while ((nextLine = reader.readLine()) != null) {
            nextLine = nextLine.trim();

            if (nextLine.startsWith("//") || nextLine.length() == 0) {
                continue;
            } else if (nextLine.startsWith(SEPARATOR_KEY)) {
                prefList.add(new Preference(SEPARATOR_KEY, group));
                continue;
            } else if (nextLine.startsWith(INFO_KEY)) {
                Preference preference = new Preference(INFO_KEY, nextLine.substring(INFO_KEY.length()).trim(), group);
                prefList.add(preference);   // "Blank" preference
                continue;
            } else if (nextLine.startsWith("##")) {

                group = null;  // End previous group
                if (nextLine.length() > 2) {
                    group = nextLine.substring(2);  // New group
                }
                continue;
            } else if (nextLine.startsWith("#")) {

                // New tab
                String[] tokens = Globals.tabPattern.split(nextLine);
                String tabLabel = tokens[0].substring(1);
                String category = tokens.length > 1 ? tokens[1] : null;

                prefList = new ArrayList<>();
                PreferenceGroup preferenceGroup = new PreferenceGroup(tabLabel, category, prefList);
                groupList.add(preferenceGroup);

                group = null;

                continue;
            } else {

                String[] tokens = Globals.tabPattern.split(nextLine);

                if (tokens[0].equals("SAM.SHOW_JUNCTION_FLANKINGREGIONS")) {
                    System.out.println();
                }

                if (tokens.length < 4) {
                    if (tokens.length == 2) {
                        // Hidden preference (not shown in editor)
                        tokens = new String[]{tokens[0], "", "", tokens[1]};
                        prefList.add(new Preference(tokens, group));
                    }

                } else {
                    prefList.add(new Preference(tokens, group));
                }
            }

        }

        return groupList;
    }

    private synchronized void storePreferences() {

        FileWriter fileWriter = null;
        try {
            fileWriter = new FileWriter(prefFile);
            PrintWriter pw = new PrintWriter(new BufferedWriter(fileWriter));

            for (Map.Entry<String, IGVPreferences> entry : preferencesMap.entrySet()) {

                if (!entry.getKey().equals(NULL_CATEGORY)) {
                    pw.println();
                    pw.println("##" + entry.getKey());
                }

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

    @Override
    public void receiveEvent(Object event) {
        if (event instanceof PreferencesChangeEvent) {
            storePreferences();
        }
    }


    static class Preference {

        String group;
        String[] tokens;

        Preference(String[] tokens, String group) {
            this.tokens = tokens;
            this.group = group;
        }

        Preference(String key, String group) {
            this(new String[]{key, null, null, null}, group);
        }

        Preference(String key, String label, String group) {
            this(new String[]{key, label, null, null}, group);
        }


        String getKey() {
            return tokens[0];
        }

        String getLabel() {
            return tokens[1];
        }

        String getType() {
            return tokens[2];
        }

        String getDefaultValue() {
            return tokens[3];
        }

        String getComment() {
            return tokens.length > 4 ? tokens[4] : null;
        }

        String getGroup() {
            return group;
        }

        String printString() {
            String str = getKey() + "\t" + getLabel() + "\t" + getType() + "\t" + getDefaultValue();
            if (getComment() != null) str += "\t" + getComment();
            return str;
        }

    }

    static class PreferenceGroup {

        String tabLabel;
        String category;
        List<Preference> preferences;

        public PreferenceGroup(String tabLabel, String category, List<Preference> preferences) {
            this.tabLabel = tabLabel;
            this.category = category;
            this.preferences = preferences;
        }
    }

//    static {
//
//        genericDefaults.put(MUTATION_INDEL_COLOR_KEY, "0,200,0");
//        genericDefaults.put(MUTATION_MISSENSE_COLOR_KEY, "170,20,240");
//        genericDefaults.put(MUTATION_NONSENSE_COLOR_KEY, "50,30,75");
//        genericDefaults.put(MUTATION_SPLICE_SITE_COLOR_KEY, "150,0,150");
//        genericDefaults.put(MUTATION_SYNONYMOUS_COLOR_KEY, "200,170,200");
//        genericDefaults.put(MUTATION_TARGETED_REGION_COLOR_KEY, "236,155,43");
//        genericDefaults.put(MUTATION_UNKNOWN_COLOR_KEY, "0,180,225");
//        //     * Nico's labels:   Truncating, Non-coding_Transcript, Other_AA_changing, Other_likely_neutral.
//        genericDefaults.put("MUTATION_Truncating_COLOR", "150,0,0");
//        genericDefaults.put("MUTATION_Non-coding_Transcript_COLOR", "0,0,150");
//        genericDefaults.put("MUTATION_Other_AA_changing_COLOR", "0,150,150");
//        genericDefaults.put("MUTATION_Other_likely_neutral_COLOR", "225,180,225");
//
//
//        genericDefaults.put(PROBE_MAPPING_KEY, "false");
//        genericDefaults.put(PROBE_MAPPING_FILE, null);
//        genericDefaults.put(USE_PROBE_MAPPING_FILE, "false");
//
//        genericDefaults.put(SHOW_REGION_BARS, "false");
//
//        genericDefaults.put(OVERLAY_MUTATION_TRACKS, "true");
//        genericDefaults.put(SHOW_ORPHANED_MUTATIONS, "true");
//        genericDefaults.put(COLOR_MUTATIONS, "false");
//        genericDefaults.put(OVERLAY_MUTATIONS_WHOLE_GENOME, "true");
//        genericDefaults.put(SHOW_SINGLE_TRACK_PANE_KEY, "false");
//        genericDefaults.put(PORT_ENABLED, "true");
//        genericDefaults.put(EXPAND_FEAUTRE_TRACKS, "false");
//        genericDefaults.put(SHOW_ATTRIBUTE_VIEWS_KEY, "true");
//        genericDefaults.put(SHOW_SINGLE_TRACK_PANE_KEY, "false");
//        genericDefaults.put(SHOW_DEFAULT_TRACK_ATTRIBUTES, "false");
//
//        genericDefaults.put(CHART_DRAW_TOP_BORDER, "false");
//        genericDefaults.put(CHART_DRAW_BOTTOM_BORDER, "false");
//        genericDefaults.put(CHART_COLOR_BORDERS, "true");
//        genericDefaults.put(CHART_DRAW_TRACK_NAME, "false");
//        genericDefaults.put(CHART_DRAW_Y_AXIS, "false");
//        genericDefaults.put(CHART_AUTOSCALE, "false");
//        genericDefaults.put(CHART_SHOW_DATA_RANGE, "true");
//        genericDefaults.put(CHART_COLOR_TRACK_NAME, "true");
//        genericDefaults.put(CHART_TRACK_HEIGHT_KEY, "40");
//        genericDefaults.put(CHART_SHOW_ALL_HEATMAP, "false");
//
//        genericDefaults.put(SAM_SHOW_DUPLICATES, "false");
//        genericDefaults.put(SAM_FILTER_DUPLICATES, "true");
//        genericDefaults.put(SAM_QUICK_CONSENSUS_MODE, "false");
//        genericDefaults.put(SAM_SHOW_SOFT_CLIPPED, "false");
//        genericDefaults.put(SAM_FLAG_UNMAPPED_PAIR, "false");
//        genericDefaults.put(SAM_AUTO_SORT, "false");
//        genericDefaults.put(SAM_SHADE_CENTER, "true");
//        genericDefaults.put(SAM_SHOW_REF_SEQ, "false");
//        genericDefaults.put(SAM_SHOW_CENTER_LINE, "true");
//        genericDefaults.put(SAM_SHOW_COV_TRACK, "true");
//        genericDefaults.put(SAM_SHADE_BASES, AlignmentTrack.ShadeBasesOption.QUALITY.toString());
//        genericDefaults.put(SAM_FILTER_ALIGNMENTS, "false");
//        genericDefaults.put(SAM_FILTER_SECONDARY_ALIGNMENTS, "false");
//        genericDefaults.put(SAM_FILTER_SUPPLEMENTARY_ALIGNMENTS, "false");
//        genericDefaults.put(SAM_FILTER_FAILED_READS, "true");
//        genericDefaults.put(SAM_DOWNSAMPLE_READS, "true");
//        genericDefaults.put(SAM_SAMPLING_WINDOW, "50");
//        genericDefaults.put(SAM_SAMPLING_COUNT, "100");
//        genericDefaults.put(SAM_BASE_QUALITY_MIN, "5");
//        genericDefaults.put(SAM_BASE_QUALITY_MAX, "20");
//        genericDefaults.put(SAM_FILTER_URL, null);
//        genericDefaults.put(SAM_HIDDEN_TAGS, "SA,MD,XA,RG");
//        genericDefaults.put(SAM_QUALITY_THRESHOLD, "0");
//        genericDefaults.put(SAM_ALLELE_THRESHOLD, "0.2f");
//        genericDefaults.put(SAM_ALLELE_USE_QUALITY, "true");
//        genericDefaults.put(SAM_MIN_INSERT_SIZE_THRESHOLD, "50");
//        genericDefaults.put(SAM_MAX_INSERT_SIZE_THRESHOLD, "1000");
//        genericDefaults.put(SAM_MIN_INSERT_SIZE_PERCENTILE, "0.5");
//        genericDefaults.put(SAM_MAX_INSERT_SIZE_PERCENTILE, "99.5");
//        genericDefaults.put(SAM_MAX_VISIBLE_RANGE, "30");
//        genericDefaults.put(SAM_COLOR_BY, "UNEXPECTED_PAIR");
//        genericDefaults.put(SAM_COLOR_BY_TAG, "");
//        genericDefaults.put(SAM_GROUP_BY_TAG, "");
//        genericDefaults.put(SAM_GROUP_BY_POS, "");
//        genericDefaults.put(SAM_SORT_BY_TAG, "");
//        genericDefaults.put(SAM_BISULFITE_CONTEXT, "CG");
//        genericDefaults.put(SAM_COMPUTE_ISIZES, "true");
//        genericDefaults.put(SAM_FLAG_ZERO_QUALITY, "true");
//        genericDefaults.put(SAM_SHOW_JUNCTION_TRACK, "false");
//        genericDefaults.put(SAM_JUNCTION_MIN_FLANKING_WIDTH, "0");
//        genericDefaults.put(SAM_JUNCTION_MIN_COVERAGE, "1");
//        genericDefaults.put(SAM_SHOW_JUNCTION_FLANKINGREGIONS, "true");
//        genericDefaults.put(SAM_NOMESEQ_ENABLED, "false");
//        genericDefaults.put(SAM_COUNT_DELETED_BASES_COVERED, "false");
//        genericDefaults.put(SAM_FLAG_LARGE_INDELS, "true");
//        genericDefaults.put(SAM_LARGE_INDELS_THRESHOLD, "1");
//        genericDefaults.put(SAM_FLAG_CLIPPING, "false");
//        genericDefaults.put(SAM_CLIPPING_THRESHOLD, "0");
//        genericDefaults.put(SAM_SORT_OPTION, "NUCLEOTIDE");
//        genericDefaults.put(SAM_GROUP_OPTION, "NONE");
//        genericDefaults.put(SAM_SHOW_GROUP_SEPARATOR, "true");
//        genericDefaults.put(SAM_COMPLETE_READS_ONLY, "false");
//        genericDefaults.put(SAM_SHOW_ALL_BASES, "false");
//
//        genericDefaults.put(SAM_REDUCED_MEMORY_MODE, "false");
//
//        genericDefaults.put(SAM_HIDE_SMALL_INDEL, "false");
//        genericDefaults.put(SAM_SMALL_INDEL_BP_THRESHOLD, "0");
//        genericDefaults.put(SAM_SHOW_INSERTION_MARKERS, "false");
//
//        genericDefaults.put(SAM_LINK_READS, "false");
//        genericDefaults.put(SAM_LINK_TAG, "READNAME");
//
//        genericDefaults.put(SAM_SHOW_ALIGNMENT_TRACK, "true");
//
//        genericDefaults.put(BYPASS_FILE_AUTO_DISCOVERY, "false");
//
//        genericDefaults.put(NORMALIZE_COVERAGE, "false");
//
//        genericDefaults.put(SHOW_GENOME_SERVER_WARNING, "true");
//
//        genericDefaults.put(SEARCH_ZOOM, "true");
//
//
//        genericDefaults.put(GENOMES_SERVER_URL, Globals.DEFAULT_GENOME_URL);
//        genericDefaults.put(OVERLAY_ATTRIBUTE_KEY, "LINKING_ID");
//        genericDefaults.put(DEFAULT_GENOME, Globals.DEFAULT_GENOME);
//
//        genericDefaults.put(USE_PROXY, "false");
//        genericDefaults.put(PROXY_AUTHENTICATE, "false");
//        genericDefaults.put(PORT_NUMBER, "60151");
//        genericDefaults.put(TRACK_HEIGHT_KEY, "15");
//        genericDefaults.put(FLANKING_REGION, "2000");
//
//        genericDefaults.put(SHOW_SEQUENCE_TRANSLATION, "false");
//        genericDefaults.put(MAX_SEQUENCE_RESOLUTION, "2");
//
//        genericDefaults.put(AUTO_UPDATE_GENOMES, "true");
//
//        genericDefaults.put(GWAS_TRACK_HEIGHT, "200");
//        genericDefaults.put(GWAS_DESCRIPTION_CACHE_SIZE, "10000");
//        genericDefaults.put(GWAS_MIN_POINT_SIZE, "3");
//        genericDefaults.put(GWAS_MAX_POINT_SIZE, "7");
//        genericDefaults.put(GWAS_USE_CHR_COLORS, "true");
//        genericDefaults.put(GWAS_SINGLE_COLOR, "false");
//        genericDefaults.put(GWAS_ALTERNATING_COLORS, "false");
//        genericDefaults.put(GWAS_PRIMARY_COLOR, "69,101,183");
//        genericDefaults.put(GWAS_SECONDARY_COLOR, "250,169,10");
//        genericDefaults.put(GWAS_SHOW_AXIS, "true");
//
//        genericDefaults.put(DEFAULT_FONT_SIZE, "10");
//        genericDefaults.put(DEFAULT_FONT_FAMILY, "Arial");
//        genericDefaults.put(DEFAULT_FONT_ATTRIBUTE, String.valueOf(Font.PLAIN));
//        genericDefaults.put(SCALE_FONTS, "false");
//
//        boolean isMac = System.getProperty("os.name").toLowerCase().startsWith("mac");
//        genericDefaults.put(ENABLE_ANTIALISING, String.valueOf(isMac));
//
//        genericDefaults.put(NAME_PANEL_WIDTH, "160");
//        genericDefaults.put(BACKGROUND_COLOR, "250,250,250");
//
//        genericDefaults.put(GENOME_SPACE_ENABLE, "true");
//        genericDefaults.put(GENOME_SPACE_DM_SERVER, "https://dm.genomespace.org/datamanager/v1.0/");
//        genericDefaults.put(GENOME_SPACE_ATM_SERVER, "https://atm.genomespace.org/atm/v1.0/");
//        genericDefaults.put(GENOME_SPACE_IDENTITY_SERVER, "https://identitydev.genomespace.org:8444/identityServer/basic");
//
//        genericDefaults.put(DB_ENABLED, "false");
//        genericDefaults.put(DB_HOST, "");
//        genericDefaults.put(DB_NAME, "");
//        genericDefaults.put(DB_PORT, "-1");
//
//
//        String defaultDataURL = Globals.DEFAULT_DATA_URL;
//        Properties properties = new Properties();
//        try {
//            properties.load(AboutDialog.class.getResourceAsStream("/resources/about.properties"));
//            String tmp = properties.getProperty("master-resource-url");
//            if (tmp != null && !tmp.startsWith("@")) {
//                defaultDataURL = tmp;
//            }
//        } catch (IOException e) {
//            log.error("Error reading dataURL property", e);
//        }
//
//        genericDefaults.put(DATA_SERVER_URL_KEY, defaultDataURL);
//
//        genericDefaults.put(CBIO_MUTATION_THRESHOLD, "1");
//        genericDefaults.put(CBIO_AMPLIFICATION_THRESHOLD, "0.9");
//        genericDefaults.put(CBIO_DELETION_THRESHOLD, "0.9");
//        genericDefaults.put(CBIO_EXPRESSION_UP_THRESHOLD, "1.0");
//        genericDefaults.put(CBIO_EXPRESSION_DOWN_THRESHOLD, "1.0");
//
//        genericDefaults.put(TOOLTIP_INITIAL_DELAY, "50");
//        genericDefaults.put(TOOLTIP_RESHOW_DELAY, "50");
//        genericDefaults.put(TOOLTIP_DISMISS_DELAY, "60000");
//        genericDefaults.put(DETAILS_BEHAVIOR_KEY, IGVCommandBar.SHOW_DETAILS_BEHAVIOR.HOVER.name());
//
//        genericDefaults.put(SHOW_SIZE_WARNING, "true");
//
//        genericDefaults.put(SKIP_VERSION, "");
//
//        genericDefaults.put(COLOR_A, "0,150,0");
//        genericDefaults.put(COLOR_C, "0,0,255");
//        genericDefaults.put(COLOR_T, "255,0,0");
//        genericDefaults.put(COLOR_G, "209,113,5");
//        genericDefaults.put(COLOR_N, ColorUtilities.colorToString(Color.gray));
//        genericDefaults.put(SAM_COLOR_A, "0,255,0");
//        genericDefaults.put(SAM_COLOR_C, "0,0,255");
//        genericDefaults.put(SAM_COLOR_T, "255,0,0");
//        genericDefaults.put(SAM_COLOR_G, "209,113,5");
//        genericDefaults.put(SAM_COLOR_N, ColorUtilities.colorToString(Color.gray.brighter()));
//
//        genericDefaults.put(HOMREF_COLOR, "235,235,235");
//        genericDefaults.put(HETVAR_COLOR, "0,0,255");
//        genericDefaults.put(HOMVAR_COLOR, "0,245,255");
//        genericDefaults.put(NOCALL_COLOR, "255,255,255");
//        genericDefaults.put(AF_REF_COLOR, "0,0,220");
//        genericDefaults.put(AF_VAR_COLOR, "255,0,0");
//
//        genericDefaults.put(VARIANT_COLOR_BY_ALLELE_FREQ, "true");
//
//        genericDefaults.put(SASHIMI_SHOW_COVERAGE, "true");
//
//        genericDefaults.put(ENABLE_GOOGLE_MENU, "false");
//        genericDefaults.put(SAVE_GOOGLE_CREDENTIALS, "true");
//
//        genericDefaults.put(DEFAULT_VISIBILITY_WINDOW, "-1");
//
//        genericDefaults.put(BLAT_URL, "http://genome.cse.ucsc.edu/cgi-bin/hgBlat");
//
//        genericDefaults.put(GENE_LIST_BED_FORMAT, "false");
//
//        genericDefaults.put(SESSION_RELATIVE_PATH, "false");
//
//        genericDefaults.put(SHOW_LOS, "true");
//
//
//        thirdGenDefaults.put(SAM_QUICK_CONSENSUS_MODE, "true");
//        thirdGenDefaults.put(SAM_DOWNSAMPLE_READS, "false");
//        thirdGenDefaults.put(SAM_MAX_VISIBLE_RANGE, "1000");
//        thirdGenDefaults.put(SAM_FLAG_LARGE_INDELS, "true");
//        thirdGenDefaults.put(SAM_LARGE_INDELS_THRESHOLD, "1");
//        thirdGenDefaults.put(SAM_FLAG_CLIPPING, "true");
//        thirdGenDefaults.put(SAM_CLIPPING_THRESHOLD, "0");
//        thirdGenDefaults.put(SAM_HIDE_SMALL_INDEL, "true");
//        thirdGenDefaults.put(SAM_SMALL_INDEL_BP_THRESHOLD, "3");
//        thirdGenDefaults.put(SAM_SHOW_INSERTION_MARKERS, "true");
//
//
//        rnaDefaults.put(SAM_MAX_VISIBLE_RANGE, "300");
//        rnaDefaults.put(SAM_SHOW_JUNCTION_TRACK, "true");
//
//
//    }
}
