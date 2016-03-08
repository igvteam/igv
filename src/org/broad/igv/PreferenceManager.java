/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv;


import org.apache.log4j.Logger;
import org.broad.igv.feature.genome.GenomeListItem;
import org.broad.igv.renderer.ColorScaleFactory;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.sam.AlignmentTrack.ShadeBasesOption;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.*;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.color.PaletteColorTable;
import org.broad.igv.ui.event.AlignmentTrackEvent;
import org.broad.igv.ui.util.PropertyManager;
import org.broad.igv.util.HttpUtils;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Manages user preferences.
 */
public class PreferenceManager implements PropertyManager {

    private static Logger log = Logger.getLogger(PreferenceManager.class);

    public static final String DEFAULT_GENOME = "hg19";

    public static final String SKIP_VERSION = "SKIP_VERSION";

    public static final String SHOW_LOS = "showLOS";

    public static final String INITIAL_TRACK_HEIGHT = "15";

    public static final String TOOLTIP_INITIAL_DELAY = "TOOLTIP.INITIAL_DELAY";
    public static final String TOOLTIP_RESHOW_DELAY = "TOOLTIP.RESHOW_DELAY";
    public static final String TOOLTIP_DISMISS_DELAY = "TOOLTIP.DISMISS_DELAY";


    public static final String CHART_DRAW_TOP_BORDER = "CHART.DRAW_TOP_BORDER";
    public static final String CHART_DRAW_BOTTOM_BORDER = "CHART.DRAW_BOTTOM_BORDER";
    public static final String CHART_COLOR_BORDERS = "CHART.COLOR_BORDERS";
    public static final String CHART_DRAW_Y_AXIS = "CHART.DRAW_AXIS";
    public static final String CHART_DRAW_TRACK_NAME = "CHART.DRAW_TRACK_NAME";
    public static final String CHART_COLOR_TRACK_NAME = "CHART.COLOR_TRACK_NAME";
    public static final String CHART_AUTOSCALE = "CHART.AUTOSCALE";
    public static final String CHART_SHOW_DATA_RANGE = "CHART.SHOW_DATA_RANGE";

    public static final String UNLOAD_ON_GENOME_CHANGE = "UNLOAD_ON_GENOME_CHANGE";

    /**
     * Added by Chantal Roth, June 25th 2012
     */
    public static final String IONTORRENT_FLOWDIST_HIDE_FIRST_HP = "IONTORRENT.FLOWDIST_HIDE_FIRST_HP";
    public static final String IONTORRENT_FLOWDIST_BINSIZE = "IONTORRENT.FLOWDIST_BINSIZE";
    public static final String IONTORRENT_FLOWDIST_CHARTTYPE = "IONTORRENT.FLOWDIST_CHARTTYPE";
    public static final String IONTORRENT_SERVER = "IONTORRENT.SERVER";
    public static final String IONTORRENT_RESULTS = "IONTORRENT.RESULTS";

    public static final String SAM_ALLELE_THRESHOLD = "SAM.ALLELE_THRESHOLD";
    public static final String SAM_ALLELE_USE_QUALITY = "SAM.ALLELE_USE_QUALITY";
    public static final String SAM_QUALITY_THRESHOLD = "SAM.QUALITY_THRESHOLD";
    public static final String SAM_MAX_INSERT_SIZE_THRESHOLD = "SAM.INSERT_SIZE_THRESHOLD";
    public static final String SAM_MIN_INSERT_SIZE_THRESHOLD = "SAM.MIN_INSERT_SIZE_THRESHOLD";
    public static final String SAM_MAX_INSERT_SIZE_PERCENTILE = "SAM.ISIZE_MAX_PERCENTILE";
    public static final String SAM_MIN_INSERT_SIZE_PERCENTILE = "SAM.MIN_ISIZE_MIN_PERCENTILE";
    public static final String SAM_AUTO_SORT = "SAM.AUTOSORT";
    public static final String SAM_SHADE_CENTER = "SAM.SHADE_CENTER";
    public static final String SAM_SHOW_CENTER_LINE = "SAM.SHOW_CENTER_LINE";
    public static final String SAM_SHOW_REF_SEQ = "SAM.SHOW_REF_SEQ";
    public static final String SAM_SHOW_COV_TRACK = "SAM.SHOW_COV_TRACK";
    public static final String SAM_SHADE_BASES = "SAM.SHADE_BASE_QUALITY";
    public static final String SAM_BASE_QUALITY_MIN = "SAM.BASE_QUALITY_MIN";
    public static final String SAM_BASE_QUALITY_MAX = "SAM.BASE_QUALITY_MAX";
    public static final String SAM_FILTER_ALIGNMENTS = "SAM.FILTER_ALIGNMENTS";
    public static final String SAM_FILTER_SECONDARY_ALIGNMENTS = "SAM.FILTER_SECONDARY_ALIGNMENTS";
    public static final String SAM_FILTER_SUPPLEMENTARY_ALIGNMENTS = "SAM.FILTER_SUPPLEMENTARY_ALIGNMENTS";
    public static final String SAM_FILTER_URL = "SAM.FILTER_URL";
    public static final String SAM_MAX_VISIBLE_RANGE = "SAM.MAX_VISIBLE_RANGE";
    public static final String SAM_SHOW_DUPLICATES = "SAM.SHOW_DUPLICATES";
    public static final String SAM_SHOW_SOFT_CLIPPED = "SAM.SHOW_SOFT_CLIPPED";
    public static final String SAM_FLAG_UNMAPPED_PAIR = "SAM.FLAG_UNMAPPED_PAIR";
    public static final String SAM_SAMPLING_COUNT = "SAM.MAX_LEVELS"; // Sampling count
    public static final String SAM_SAMPLING_WINDOW = "SAM.SAMPLING_WINDOW";
    public static final String SAM_DOWNSAMPLE_READS = "SAM.DOWNSAMPLE_READS";
    public static final String SAM_SORT_OPTION = "SAM.SORT_OPTION";
    public static final String SAM_SHOW_ALL_BASES = "SAM.SHOW_ALL_BASES";

    public static final String SAM_COLOR_BY = "SAM.COLOR_BY";
    public static final String SAM_COLOR_BY_TAG = "SAM.COLOR_BY_TAG";
    public static final String SAM_SORT_BY_TAG = "SAM.SORT_BY_TAG";
    public static final String SAM_GROUP_BY_TAG = "SAM.GROUP_BY_TAG";
    public static final String SAM_BISULFITE_CONTEXT = "SAM.BISULFITE_CONTEXT";
    public static final String SAM_FILTER_FAILED_READS = "SAM.FILTER_FAILED_READS";
    public static final String SAM_COMPUTE_ISIZES = "SAM.COMPUTE_ISIZES";
    public static final String SAM_FLAG_ZERO_QUALITY = "SAM.FLAG_ZERO_QUALITY";
    //dhmay adding 20110208
    public static final String SAM_SHOW_JUNCTION_TRACK = "SAM.SHOW_JUNCTION_TRACK";
    public static final String SAM_JUNCTION_MIN_FLANKING_WIDTH = "SAM.JUNCTION_MIN_FLANKING_WIDTH";
    public static final String SAM_JUNCTION_MIN_COVERAGE = "SAM.JUNCTION_MIN_COVERAGE";
    //dhmay adding 20120731
    public static final String SAM_SHOW_JUNCTION_FLANKINGREGIONS = "SAM.SHOW_JUNCTION_FLANKINGREGIONS";

    public static final String SAM_NOMESEQ_ENABLED = "SAM.NOMESEQ_ENABLED";
    public static final String SAM_COUNT_DELETED_BASES_COVERED = "SAM.COUNT_DELETED_BASES_COVERED";

    public static final String SAM_FLAG_LARGE_INSERTIONS = "SAM.FLAG_LARGE_INSERTIONS";
    public static final String SAM_LARGE_INSERTIONS_THRESHOLD = "SAM.LARGE_INSERTIONS_THRESOLD";

    public static final String SAM_SHOW_GROUP_SEPARATOR = "SAM.SHOW_GROUP_SEPARATOR";
    public static final String SAM_COMPLETE_READS_ONLY = "SAM.COMPLETE_READS_ONLY";

    public static final String SAM_REDUCED_MEMORY_MODE = "SAM.REDUCED_MEMORY_MODE";
    public static final String SAM_MIN_INDEL_SIZE = "SAM.MIN_INDEL_SIZE";

    public static final String SAM_SHOW_ALIGNMENT_TRACK = "SAM.SHOW_ALIGNMENT_TRACK";

    public static final String COLOR_A = "COLOR.A";
    public static final String COLOR_C = "COLOR.C";
    public static final String COLOR_T = "COLOR.T";
    public static final String COLOR_G = "COLOR.G";
    public static final String COLOR_N = "COLOR.N";
    public static final String SAM_COLOR_A = "SAM.COLOR.A";
    public static final String SAM_COLOR_C = "SAM.COLOR.C";
    public static final String SAM_COLOR_T = "SAM.COLOR.T";
    public static final String SAM_COLOR_G = "SAM.COLOR.G";
    public static final String SAM_COLOR_N = "SAM.COLOR.N";

    public static final String HOMREF_COLOR = "HOMREF.COLOR";
    public static final String HETVAR_COLOR = "HETVAR.COLOR";
    public static final String HOMVAR_COLOR = "HOMVAR.COLOR";
    public static final String NOCALL_COLOR = "NOCALL.COLOR";
    public static final String AF_REF_COLOR = "AF_REF.COLOR";
    public static final String AF_VAR_COLOR = "AF_VAR.COLOR";

    public static final String EXPAND_FEAUTRE_TRACKS = "EXPAND_FEATURE_TRACKS";
    public static final String PORT_ENABLED = "PORT_ENABLED";
    public static final String PORT_NUMBER = "PORT_NUMBER";
    public static final String COLOR_SCALE_KEY = "COLOR_SCALE_";
    public static final String FRAME_BOUNDS_KEY = "IGV.Bounds";
    public static final String FRAME_STATE_KEY = "IGV.Frame.ExtendedState";
    public static final String RECENT_SESSION_KEY = "IGV.Session.recent.sessions";
    public static final String TRACK_HEIGHT_KEY = "IGV.track.height";
    public static final String CHART_TRACK_HEIGHT_KEY = "IGV.chart.track.height";
    public static final String CHART_SHOW_ALL_HEATMAP = "CHART.SHOW_ALL_HEATMAP";
    public static final String SHOW_MISSING_DATA_KEY = "IGV.track.show.missing.data";
    public static final String SHOW_ATTRIBUTE_VIEWS_KEY = "IGV.track.show.attribute.views";
    public static final String SHOW_SINGLE_TRACK_PANE_KEY = "IGV.single.track.pane";
    public static final String GENOMES_SERVER_URL = "IGV.genome.sequence.dir";
    public static final String JOIN_ADJACENT_SEGMENTS_KEY = "IGV.join.adjacent.segments";
    public static final String SHOW_REGION_BARS = "SHOW_REGION_BARS";
    public static final String LAST_EXPORTED_REGION_DIRECTORY = "LAST_EXPORTED_REGION_DIRECTORY";
    final static public String LAST_TRACK_DIRECTORY = "LAST_TRACK_DIRECTORY";
    final static public String LAST_SNAPSHOT_DIRECTORY = "LAST_SNAPSHOT_DIRECTORY";
    final static public String LAST_GENOME_IMPORT_DIRECTORY = "LAST_GENOME_IMPORT_DIRECTORY";
    final static public String LAST_SESSION_DIRECTORY = "LAST_SESSION_DIRECTORY";
    final static public String DEFAULT_GENOME_KEY = "DEFAULT_GENOME_KEY";
    final static public String LAST_CHROMOSOME_VIEWED_KEY = "LAST_CHROMOSOME_VIEWED_KEY";
    final static public String HISTORY_DELIMITER = ";";
    final static public String GENOME_ID_DISPLAY_LIST_KEY = "GENOME_LIST";
    final static public String DETAILS_BEHAVIOR_KEY = "DETAILS_BEHAVIOR";
    final static public String DEFAULT_VISIBILITY_WINDOW = "DEFAULT_VISIBILITY_WINDOW";

    public static final String MUTATION_COLOR_TABLE = "MUTATION_COLOR_TABLE";
    public static final String MUTATION_INDEL_COLOR_KEY = "MUTATION_INDEL_COLOR_KEY";
    public static final String MUTATION_MISSENSE_COLOR_KEY = "MUTATION_MISSENSE_COLOR_KEY";
    public static final String MUTATION_NONSENSE_COLOR_KEY = "MUTATION_NONSENSE_COLOR_KEY";
    public static final String MUTATION_SPLICE_SITE_COLOR_KEY = "MUTATION_SPLICE_SITE_COLOR_KEY";
    public static final String MUTATION_SYNONYMOUS_COLOR_KEY = "MUTATION_SYNONYMOUS_COLOR_KEY";
    public static final String MUTATION_TARGETED_REGION_COLOR_KEY = "MUTATION_TARGETED_REGION_COLOR_KEY";
    public static final String MUTATION_UNKNOWN_COLOR_KEY = "MUTATION_UNKNOWN_COLOR_KEY";
    public static final String OVERLAY_MUTATION_TRACKS = "OVERLAY_TRACKS_KEY";
    public static final String SHOW_ORPHANED_MUTATIONS = "SHOW_ORPHANED_MUTATIONS";
    public static final String OVERLAY_ATTRIBUTE_KEY = "OVERLAY_ATTRIBUTE_KEY";
    public static final String OVERLAY_MUTATIONS_WHOLE_GENOME = "OVERLAY_MUTATIONS_WHOLE_GENOME";
    public static final String COLOR_MUTATIONS = "COVER_OVERLAY_KEY";
    public static final String TRACK_ATTRIBUTE_NAME_KEY = "TRACK_ATTRIBUTE_NAME_KEY";
    public static final String DATA_SERVER_URL_KEY = "MASTER_RESOURCE_FILE_KEY";
    //public static final String CHECKED_RESOURCES_KEY = "CHECKED_RESOURCES_KEY";
    public static final String DEFINE_GENOME_INPUT_DIRECTORY_KEY = "DEFINE_GENOME_INPUT_DIRECTORY_KEY";

    public static final String VARIANT_COLOR_BY_ALLELE_FREQ = "VARIANT_COLOR_BY_ALLELE_FREQ";

    public static final String PROBE_MAPPING_KEY = "PROBE_MAPPING_KEY";
    public static final String PROBE_MAPPING_FILE = "PROBE_MAPPING_FILE";
    public static final String USE_PROBE_MAPPING_FILE = "USE_PROBE_MAPPING_FILE";

    public static final String SEARCH_ZOOM = "SEARCH_ZOOM";
    public static final String BYPASS_FILE_AUTO_DISCOVERY = "BYPASS_FILE_AUTO_DISCOVERY";
    public static final String NORMALIZE_COVERAGE = "NORMALIZE_COVERAGE";
    public static final String SHOW_EXPAND_ICON = "SHOW_EXPAND_ICON";
    public static final String SHOW_DEFAULT_TRACK_ATTRIBUTES = "SHOW_DEFAULT_TRACK_ATTRIBUTES";

    public static final String SHOW_SIZE_WARNING = "SHOW_SIZE_WARNING";
    public static final String SHOW_GENOME_SERVER_WARNING = "SHOW_GENOME_SERVER_WARNING";

    public static final String USE_PROXY = "PROXY.USE";
    public static final String PROXY_HOST = "PROXY.HOST";
    public static final String PROXY_PORT = "PROXY.PORT";
    public static final String PROXY_AUTHENTICATE = "PROXY.AUTHENTICATE";
    public static final String PROXY_NTLM = "PROXY.NTLM";
    public static final String PROXY_USER = "PROXY.USERNAME";
    public static final String PROXY_PW = "PROXY.PW";
    public static final String PROXY_TYPE = "PROXY.TYPE";
    public static final String PROXY_WHITELIST = "PROXY.WHITELIST";

    public static final String KNOWN_SNPS = "KNOWN_SNPS_FILE";

    public static final String FLANKING_REGION = "FLANKING_REGION";

    public static final String SHOW_SEQUENCE_TRANSLATION = "SHOW_SEQUENCE_TRANSLATION";
    public static final String MAX_SEQUENCE_RESOLUTION = "MAX_SEQUENCE_RESOLUTION";

    public static final String AUTO_UPDATE_GENOMES = "AUTO_UPDATE_GENOMES";

    public static final String GWAS_TRACK_HEIGHT = "GWAS_TRACK_HEIGHT";
    public static final String GWAS_DESCRIPTION_CACHE_SIZE = "GWAS_DESCRIPTION_CACHE_SIZE";
    public static final String GWAS_MIN_POINT_SIZE = "GWAS_MIN_POINT_SIZE";
    public static final String GWAS_MAX_POINT_SIZE = "GWAS_MAX_POINT_SIZE";
    public static final String GWAS_USE_CHR_COLORS = "GWAS_USE_CHR_COLORS";
    public static final String GWAS_SINGLE_COLOR = "GWAS_SINGLE_COLOR";
    public static final String GWAS_ALTERNATING_COLORS = "GWAS_ALTERNATING_COLORS";
    public static final String GWAS_PRIMARY_COLOR = "GWAS_PRIMARY_COLOR";
    public static final String GWAS_SECONDARY_COLOR = "GWAS_SECONDARY_COLOR";
    public static final String GWAS_SHOW_AXIS = "GWAS_SHOW_AXIS";

    public static final String DEFAULT_FONT_SIZE = "DEFAULT_FONT_SIZE";
    public static final String DEFAULT_FONT_FAMILY = "DEFAULT_FONT_FAMILY";
    public static final String DEFAULT_FONT_ATTRIBUTE = "DEFAULT_FONT_ATTRIBUTE";
    public static final String ENABLE_ANTIALISING = "ENABLE_ANTIALIASING";
    public static final String SCALE_FONTS = "SCALE_FONTS";

    public static final String NAME_PANEL_WIDTH = "NAME_PANEL_WIDTH";
    public static final String BACKGROUND_COLOR = "BACKGROUND_COLOR";

    public static final String GENOME_SPACE_ENABLE = "GENOME_SPACE_ENABLE";
    public static final String GENOME_SPACE_DM_SERVER = "GENOME_SPACE_DM_SERVER";
    public static final String GENOME_SPACE_ATM_SERVER = "GENOME_SPACE_ATM_SERVER";
    public static final String GENOME_SPACE_IDENTITY_SERVER = "GENOME_SPACE_IDENTITY_SERVER";

    public static final String ENABLE_EXOME_BUTTON = "ENABLE_EXOME_BUTTON";

    public static final String CBIO_MUTATION_THRESHOLD = "CBIO_MUTATION_THRESHOLD";
    public static final String CBIO_AMPLIFICATION_THRESHOLD = "CBIO_AMPLIFICATION_THRESHOLD";
    public static final String CBIO_DELETION_THRESHOLD = "CBIO_DELETION_THRESHOLD";
    public static final String CBIO_EXPRESSION_UP_THRESHOLD = "CBIO_EXPRESSION_UP_THRESHOLD";
    public static final String CBIO_EXPRESSION_DOWN_THRESHOLD = "CBIO_EXPRESSION_DOWN_THRESHOLD";


    public static final String DB_ENABLED = "DB_ENABLED";
    public static final String DB_HOST = "DB_HOST";
    public static final String DB_NAME = "DB_NAME";
    public static final String DB_PORT = "DB_PORT";
    public static final String DEFAULT_GENOME_URL = "http://igv.broadinstitute.org/genomes/genomes.txt";
    public static final String DEFAULT_DATA_URL = "http://data.broadinstitute.org/igvdata/$$_dataServerRegistry.txt";

    public static final String IGV_PLUGIN_LIST_KEY = "IGV_PLUGIN_LIST";

    public static final String SASHIMI_SHOW_COVERAGE = "SASHIMI.SHOW_COVERAGE";

    public static final String GOOGLE_API_KEY = "GOOGLE_API_KEY";
    public static final String ENABLE_GOOGLE_MENU = "ENABLE_GOOGLE_MENU";
    public static final String SAVE_GOOGLE_CREDENTIALS = "SAVE_GOOGLE_CREDENTIALS";

    public static final String BLAT_URL = "BLAT_URL";

    public static final String GENE_LIST_BED_FORMAT = "GENE_LIST_BED_FORMAT";

    public static final String SESSION_RELATIVE_PATH = "SESSION.RELATIVE_PATH";

    IGVPreferences preferences;
    Map<String, String> defaultValues;


    /**
     * Cache of preference values.  Profiling reveals that Preferences.get()
     * is taking huge amounts of time.  There are hundreds of thousands of
     * calls to this to get the track height,  this is possibly a bad design
     * decision, however caching the preference values solves the performance
     * problem for now.
     */
    private Map<String, Boolean> booleanCache = new Hashtable();
    private Map<String, Object> objectCache = new Hashtable();
    private Map<TrackType, ContinuousColorScale> colorScaleCache = new Hashtable();

    private PaletteColorTable mutationColorScheme = null;


    public static PreferenceManager getInstance() {

        return instance;
    }

    private static PreferenceManager instance = new PreferenceManager();

    private PreferenceManager() {
        preferences = new IGVPreferences();
        initDefaultValues();
    }


    public String get(String key, String defaultString) {
        key = key.trim();
        return preferences.get(key, defaultString);
    }

    public String get(String key) {
        key = key.trim();
        return get(key, defaultValues.get(key));
    }

    public boolean hasExplicitValue(String key) {
        key = key.trim();
        return preferences.userPreferences.containsKey(key);
    }

    /**
     * Get the default value for the specified key.
     * May be null.
     *
     * @param key
     * @return
     */
    public String getDefaultValue(String key) {
        key = key.trim();
        return defaultValues.get(key);
    }


    /**
     * Return the preference as a boolean value.
     *
     * @param key
     * @return
     */
    public boolean getAsBoolean(String key) {
        key = key.trim();
        Boolean boolValue = booleanCache.get(key);
        if (boolValue == null) {
            String value = get(key);
            if (value == null) {
                log.error("No default value for: " + key);
                return false;
            }
            boolValue = new Boolean(get(key, value));
            booleanCache.put(key, boolValue);
        }
        return boolValue.booleanValue();
    }

    /**
     * Return the preference as an integer.
     *
     * @param key
     * @return
     */
    public int getAsInt(String key) {
        key = key.trim();
        Number value = (Number) objectCache.get(key);
        if (value == null) {
            String defValue = get(key);
            if (defValue == null) {
                log.error("No default value for: " + key);
                return 0;
            }
            value = new Integer(get(key, defValue));
            objectCache.put(key, value);
        }
        return value.intValue();
    }

    /**
     * Return the preference as a color.
     *
     * @param key
     * @return
     */
    public Color getAsColor(String key) {
        key = key.trim();
        Color value = (Color) objectCache.get(key);
        if (value == null) {
            String defValue = get(key);
            if (defValue == null) {
                log.error("No default value for: " + key);
                return Color.white;
            }
            value = ColorUtilities.stringToColor(defValue);
            objectCache.put(key, value);
        }
        return value;
    }

    /**
     * Return the preference as an float.
     *
     * @param key
     * @return
     */
    public float getAsFloat(String key) {
        key = key.trim();
        Number value = (Number) objectCache.get(key);
        if (value == null) {
            String defValue = get(key);
            if (defValue == null) {
                log.error("No default value for: " + key);
                return 0;
            }
            value = new Float(get(key, defValue));
            objectCache.put(key, value);
        }
        return value.floatValue();
    }


    public void mergePreferences(Map<String, String> newPrefs) {
        for (Map.Entry<String, String> entry : newPrefs.entrySet()) {
            String key = entry.getKey();
            if (!newPrefs.containsKey(key)) {
                put(key, entry.getValue());
            }
        }
    }


    /**
     * Update any cached values with the new key/value pair
     *
     * @param key
     * @param value
     */
    private void updateCaches(String key, String value) {
        key = key.trim();
        if (booleanCache.containsKey(key)) {
            booleanCache.put(key, new Boolean(value));
        }
        colorScaleCache.remove(key);
        objectCache.remove(key);
    }

    private void clearCaches() {
        colorScaleCache.clear();
        booleanCache.clear();
        objectCache.clear();
    }

    public void put(String key, String value) {
        key = key.trim();
        preferences.put(key, value);
        updateCaches(key, value);
    }

    public void put(String key, boolean b) {
        key = key.trim();
        String value = String.valueOf(b);
        preferences.put(key, value);
        updateCaches(key, value);
    }


    public void putAll(Map<String, String> updatedPrefs) {


        for (Map.Entry<String, String> entry : updatedPrefs.entrySet()) {
            if (entry.getValue() == null || entry.getValue().trim().length() == 0) {
                remove(entry.getKey());

            } else {
                put(entry.getKey(), entry.getValue());
            }
        }

        checkForAlignmentChanges(updatedPrefs);


        clearCaches();

    }

    private void checkForAlignmentChanges(Map<String, String> updatedPreferenceMap) {

        if (IGV.hasInstance()) {

            final IGV igv = IGV.getInstance();

            boolean reloadSAM = false;
            for (String key : SAM_RELOAD_KEYS) {
                if (updatedPreferenceMap.containsKey(key)) {
                    reloadSAM = true;
                    break;
                }

            }

            if (reloadSAM) {
                if (updatedPreferenceMap.containsKey(PreferenceManager.SAM_MAX_VISIBLE_RANGE)) {
                    igv.notifyAlignmentTrackEvent(this, AlignmentTrackEvent.Type.VISIBILITY_WINDOW);
                }
                igv.notifyAlignmentTrackEvent(this, AlignmentTrackEvent.Type.RELOAD);
            }
            if (updatedPreferenceMap.containsKey(PreferenceManager.SAM_ALLELE_THRESHOLD)) {
                igv.notifyAlignmentTrackEvent(this, AlignmentTrackEvent.Type.ALLELE_THRESHOLD);
            }
        }

    }


    public void remove(String key) {
        preferences.remove(key);
        booleanCache.remove(key);
        objectCache.remove(key);
        colorScaleCache.remove(key);
    }


    public void clear() {
        preferences.clear();
        colorScaleCache.clear();
        booleanCache.clear();
        objectCache.clear();
    }


    public String getGenomeListURL() {
        return get(PreferenceManager.GENOMES_SERVER_URL);
    }

    public void overrideGenomeServerURL(String url) {
        preferences.putOverride(GENOMES_SERVER_URL, url);
    }


    /**
     * @param bounds
     */
    public void setApplicationFrameBounds(Rectangle bounds) {

        StringBuffer buffer = new StringBuffer();
        buffer.append(bounds.x);
        buffer.append(",");
        buffer.append(bounds.y);
        buffer.append(",");
        buffer.append(bounds.width);
        buffer.append(",");
        buffer.append(bounds.height);
        put(FRAME_BOUNDS_KEY, buffer.toString());
    }

    /**
     * @return
     */
    public Rectangle getApplicationFrameBounds() {

        Rectangle bounds = null;

        // Set the application's previous location and size
        String applicationBounds = preferences.get(FRAME_BOUNDS_KEY, null);

        if (applicationBounds != null) {
            String[] values = applicationBounds.split(",");
            int x = Integer.parseInt(values[0]);
            int y = Integer.parseInt(values[1]);
            int width = Integer.parseInt(values[2]);
            int height = Integer.parseInt(values[3]);
            bounds = new Rectangle(x, y, width, height);
        }
        return bounds;
    }

    /**
     * @param directory
     */
    public void setLastExportedRegionDirectory(File directory) {

        put(LAST_EXPORTED_REGION_DIRECTORY, directory.getAbsolutePath());
    }

    /**
     * @return
     */
    public File getLastExportedRegionDirectory() {

        File exportedRegionDirectory = null;

        String lastFilePath = get(LAST_EXPORTED_REGION_DIRECTORY, null);

        if (lastFilePath != null) {

            // Create the exported region directory
            exportedRegionDirectory = new File(lastFilePath);
        }

        return exportedRegionDirectory;
    }

    /**
     * @param directory
     */
    public void setLastSnapshotDirectory(File directory) {

        put(LAST_SNAPSHOT_DIRECTORY, directory.getAbsolutePath());
    }

    /**
     * @return
     */
    public File getLastSnapshotDirectory() {

        File snapshotDirectory = null;

        String lastFilePath = get(LAST_SNAPSHOT_DIRECTORY, null);

        if (lastFilePath != null) {

            // Create the snapshot directory
            snapshotDirectory = new File(lastFilePath);
        }

        return snapshotDirectory;
    }

    /**
     * @param directory
     */
    public void setDefineGenomeInputDirectory(File directory) {

        put(DEFINE_GENOME_INPUT_DIRECTORY_KEY, directory.getAbsolutePath());
    }

    /**
     * @return
     */
    public File getDefineGenomeInputDirectory() {

        File directory = null;

        String lastFilePath = get(DEFINE_GENOME_INPUT_DIRECTORY_KEY, DirectoryManager.getUserDirectory().getAbsolutePath());

        if (lastFilePath != null) {
            directory = new File(lastFilePath);
        }

        return directory;
    }

    /**
     * @param directory
     */
    public void setLastGenomeImportDirectory(File directory) {

        put(LAST_GENOME_IMPORT_DIRECTORY, directory.getAbsolutePath());
    }

    /**
     * @return
     */
    public File getLastGenomeImportDirectory() {

        File genomeImportDirectory = null;

        String lastFilePath = get(LAST_GENOME_IMPORT_DIRECTORY, DirectoryManager.getUserDirectory().getAbsolutePath());

        if (lastFilePath != null) {
            genomeImportDirectory = new File(lastFilePath);
        }

        return genomeImportDirectory;
    }


    /**
     * @param directory
     */
    public void setLastSessionDirectory(File directory) {

        put(LAST_SESSION_DIRECTORY, directory.getAbsolutePath());
    }

    /**
     * @return
     */
    public File getLastSessionDirectory() {

        File sessionDirectory = null;

        String lastFilePath = get(LAST_SESSION_DIRECTORY, null);

        if (lastFilePath != null) {

            // Create the session directory
            sessionDirectory = new File(lastFilePath);
        }

        return sessionDirectory;
    }


    /**
     * @param recentSessions
     */
    public void setRecentSessions(String recentSessions) {
        put(RECENT_SESSION_KEY, recentSessions);
    }


    public String getRecentSessions() {

        return get(RECENT_SESSION_KEY, null);
    }

    public String getDataServerURL() {
        String masterResourceFile = get(DATA_SERVER_URL_KEY);
        return masterResourceFile;
    }

    /**
     * Temporarily override the data server url with the supplied value.  This override will persist for
     * the duration of the session, or until the user explicitly changes it.
     *
     * @param url
     */
    public void overrideDataServerURL(String url) {
        preferences.putOverride(DATA_SERVER_URL_KEY, url);
        clearCaches();
    }

    /**
     * Temporarily override a preference.   This override will persist for
     * the duration of the session, or until the user explicitly changes it.
     *
     * @param key
     * @param value
     */
    public void override(String key, String value) {
        preferences.putOverride(key, value);

        Map<String, String> updatedPrefs = new HashMap<String, String>();
        updatedPrefs.put(key, value);
        checkForAlignmentChanges(updatedPrefs);


        clearCaches();
    }

    public void loadOverrides(String overridePropertyFilePath) {
        preferences.loadOverrides(overridePropertyFilePath);
        clearCaches();
    }

    public void setShowAttributeView(boolean isShowable) {
        put(PreferenceManager.SHOW_ATTRIBUTE_VIEWS_KEY, Boolean.toString(isShowable));
    }


    public void setLastChromosomeViewed(String chromosome) {
        put(LAST_CHROMOSOME_VIEWED_KEY, chromosome);
    }


    public String getLastChromosomeViewed() {
        String chromosome = get(LAST_CHROMOSOME_VIEWED_KEY, "chr1");
        return chromosome;
    }


    public void setDefaultGenome(String genomeId) {
        put(DEFAULT_GENOME_KEY, genomeId);
    }


    public String getDefaultGenome() {

        String genome = get(DEFAULT_GENOME_KEY, DEFAULT_GENOME);
        return genome;
    }

    public void setLastTrackDirectory(File directory) {
        String lastDirectory = directory.isDirectory() ? directory.getAbsolutePath() : directory.getParent();
        put(LAST_TRACK_DIRECTORY, lastDirectory);
    }

    public File getLastTrackDirectory() {

        String lastDirectoryPath = get(LAST_TRACK_DIRECTORY, null);

        File lastDirectoryFile = null;
        if (lastDirectoryPath != null) {
            lastDirectoryFile = new File(lastDirectoryPath);
        }

        return lastDirectoryFile;
    }


    /**
     * Set the color scheme for the track type.  Its unfortunate that this is a public
     * method,  color schemes are managed by ColorScaleFactory and that is the
     * only object that should call this method.
     *
     * @param type
     * @param colorScale
     */
    public void setColorScale(TrackType type, ContinuousColorScale colorScale) {
        String colorScaleString = colorScale.asString();
        put(COLOR_SCALE_KEY + type.toString(), colorScaleString);
        colorScaleCache.put(type, colorScale);
    }

    /**
     * Get the default color scheme for the track type.  Only specific TrackTypes have default schemes.  The scale
     * returned is marked as "default".
     *
     * @param type
     * @return
     */


    public ContinuousColorScale getColorScale(TrackType type) {
        if (type == null) {
            return null;
        }

        ContinuousColorScale scale = colorScaleCache.get(type);

        if (scale == null && scaledTypes.contains(type)) {
            String colorScaleString = get(COLOR_SCALE_KEY + type.toString(), null);
            if (colorScaleString != null) {
                scale = (ContinuousColorScale) ColorScaleFactory.getScaleFromString(colorScaleString);
            } else {
                scale = getDefaultColorScale(type);
            }
            if (scale != null) {
                scale.setDefault(true);
                colorScaleCache.put(type, scale);
            }
        }
        return scale;
    }


    static Set<String> scaledTypes = new HashSet(Arrays.asList(
            TrackType.LOH, TrackType.RNAI, TrackType.POOLED_RNAI, TrackType.DNA_METHYLATION,
            TrackType.GENE_EXPRESSION, TrackType.COPY_NUMBER, TrackType.ALLELE_SPECIFIC_COPY_NUMBER, TrackType.CNV));


    /**
     * Return the default color scale.  This si the scale for track type "generic",
     * as well as any track type without a specific scale.
     *
     * @param type
     * @return
     */
    public ContinuousColorScale getDefaultColorScale(TrackType type) {
        switch (type) {
            case LOH:
                return new ContinuousColorScale(0, -1, 0, 1, Color.red, UIConstants.LIGHT_YELLOW, Color.blue);
            case RNAI:
            case POOLED_RNAI:
                ContinuousColorScale cs = new ContinuousColorScale(0, -3, 0, 3, Color.red, Color.white, Color.blue);
                cs.setNoDataColor(new Color(225, 225, 225));
                return cs;

            case DNA_METHYLATION:
                cs = new ContinuousColorScale(0, 1, Color.BLUE, Color.RED);
                cs.setNoDataColor(Color.WHITE);
                return cs;

            case GENE_EXPRESSION:
                cs = new ContinuousColorScale(-0.1, -1.5, 0.1, 1.5, Color.BLUE, Color.WHITE, Color.RED);
                cs.setNoDataColor(new Color(225, 225, 225));
                return cs;

            case COPY_NUMBER:
            case ALLELE_SPECIFIC_COPY_NUMBER:
            case CNV:
                return new ContinuousColorScale(-0.1, -1.5, 0.1, 1.5, Color.BLUE, Color.WHITE, Color.RED);

            default:
                return null;
        }
    }

    /**
     * Original labels:  Indel, Missense, Nonsesne, Splice_site, Synonymous, Targetd_Region, Unknown
     * Nico's labels:   Synonymous, Missense, Truncating, Non-coding_Transcript, Other_AA_changing, Other_likely_neutral.
     */
    public void resetMutationColorScheme() {

        remove(MUTATION_INDEL_COLOR_KEY);
        remove(MUTATION_MISSENSE_COLOR_KEY);
        remove(MUTATION_NONSENSE_COLOR_KEY);
        remove(MUTATION_SPLICE_SITE_COLOR_KEY);
        remove(MUTATION_SYNONYMOUS_COLOR_KEY);
        remove(MUTATION_TARGETED_REGION_COLOR_KEY);
        remove(MUTATION_UNKNOWN_COLOR_KEY);
        remove("MUTATION_Truncating_COLOR");
        remove("MUTATION_Non-coding_Transcript_COLOR");
        remove("MUTATION_Other_AA_changing_COLOR");
        remove("MUTATION_Other_likely_neutral_COLOR");
        remove(MUTATION_COLOR_TABLE);
    }


    /**
     * Original labels:  Indel, Missense, Nonsesne, Splice_site, Synonymous, Targetd_Region, Unknown
     * Nico's labels:   Synonymous, Missense, Truncating, Non-coding_Transcript, Other_AA_changing, Other_likely_neutral.
     * Combined: Indel, Missense, Nonsesne, Splice_site, Synonymous, Targetd_Region, Unknown, Truncating,
     * Non-coding_Transcript, Other_AA_changing, Other_likely_neutral
     */
    public synchronized PaletteColorTable getMutationColorScheme() {
        if (mutationColorScheme == null) {
            String colorTableString = get(MUTATION_COLOR_TABLE);
            if (colorTableString != null) {
                PaletteColorTable pallete = new PaletteColorTable();
                pallete.restoreMapFromString(colorTableString);
                mutationColorScheme = pallete;
            } else {
                mutationColorScheme = getLegacyMutationColorScheme();
            }
        }
        return mutationColorScheme;
    }

    private PaletteColorTable getLegacyMutationColorScheme() {
        String indelColor = get(MUTATION_INDEL_COLOR_KEY);
        String missenseColor = get(MUTATION_MISSENSE_COLOR_KEY);
        String nonsenseColor = get(MUTATION_NONSENSE_COLOR_KEY);
        String spliceSiteColor = get(MUTATION_SPLICE_SITE_COLOR_KEY);
        String synonymousColor = get(MUTATION_SYNONYMOUS_COLOR_KEY);
        String targetedRegionColor = get(MUTATION_TARGETED_REGION_COLOR_KEY);
        String unknownColor = get(MUTATION_UNKNOWN_COLOR_KEY);

        PaletteColorTable colorTable = new PaletteColorTable();
        if ((indelColor != null) && (missenseColor != null) && (nonsenseColor != null) &&
                (spliceSiteColor != null) &&
                (synonymousColor != null) &&
                (targetedRegionColor != null) &&
                (unknownColor != null)) {

            Color color1 = ColorUtilities.stringToColor(indelColor);
            colorTable.put("Indel", color1);

            Color color2 = ColorUtilities.stringToColor(missenseColor);
            colorTable.put("Missense", color2);

            Color color3 = ColorUtilities.stringToColor(nonsenseColor);
            colorTable.put("Nonsense", color3);

            Color color4 = ColorUtilities.stringToColor(spliceSiteColor);
            colorTable.put("Splice_site", color4);

            Color color5 = ColorUtilities.stringToColor(synonymousColor);
            colorTable.put("Synonymous", color5);

            Color color6 = ColorUtilities.stringToColor(targetedRegionColor);
            colorTable.put("Targeted_Region", color6);

            Color color7 = ColorUtilities.stringToColor(unknownColor);
            colorTable.put("Unknown", color7);

            // Nicos extensions
            String[] nicosCats = {"Truncating", "Non-coding_Transcript", "Other_AA_changing", "Other_likely_neutral"};
            for (String cat : nicosCats) {
                String key = "MUTATION_" + cat + "_COLOR";
                colorTable.put(cat, ColorUtilities.stringToColor(get(key)));
            }
        }

        return colorTable;
    }


    /**
     * Immediately clear all proxy settings.
     */
    public void clearProxySettings() {

        remove(PreferenceManager.USE_PROXY);
        remove(PreferenceManager.PROXY_HOST);
        remove(PreferenceManager.PROXY_PORT);
        remove(PreferenceManager.PROXY_AUTHENTICATE);
        remove(PreferenceManager.PROXY_USER);
        remove(PreferenceManager.PROXY_PW);
        remove(PreferenceManager.PROXY_TYPE);
        remove(PreferenceManager.PROXY_WHITELIST);
        HttpUtils.getInstance().updateProxySettings();
    }


    private void initDefaultValues() {

        defaultValues = new HashMap();

        defaultValues.put(MUTATION_INDEL_COLOR_KEY, "0,200,0");
        defaultValues.put(MUTATION_MISSENSE_COLOR_KEY, "170,20,240");
        defaultValues.put(MUTATION_NONSENSE_COLOR_KEY, "50,30,75");
        defaultValues.put(MUTATION_SPLICE_SITE_COLOR_KEY, "150,0,150");
        defaultValues.put(MUTATION_SYNONYMOUS_COLOR_KEY, "200,170,200");
        defaultValues.put(MUTATION_TARGETED_REGION_COLOR_KEY, "236,155,43");
        defaultValues.put(MUTATION_UNKNOWN_COLOR_KEY, "0,180,225");
        //     * Nico's labels:   Truncating, Non-coding_Transcript, Other_AA_changing, Other_likely_neutral.
        defaultValues.put("MUTATION_Truncating_COLOR", "150,0,0");
        defaultValues.put("MUTATION_Non-coding_Transcript_COLOR", "0,0,150");
        defaultValues.put("MUTATION_Other_AA_changing_COLOR", "0,150,150");
        defaultValues.put("MUTATION_Other_likely_neutral_COLOR", "225,180,225");


        defaultValues.put(PROBE_MAPPING_KEY, "false");
        defaultValues.put(PROBE_MAPPING_FILE, null);
        defaultValues.put(USE_PROBE_MAPPING_FILE, "false");

        defaultValues.put(SHOW_REGION_BARS, "false");
        defaultValues.put(JOIN_ADJACENT_SEGMENTS_KEY, "false");

        defaultValues.put(OVERLAY_MUTATION_TRACKS, "true");
        defaultValues.put(SHOW_ORPHANED_MUTATIONS, "true");
        defaultValues.put(COLOR_MUTATIONS, "false");
        defaultValues.put(OVERLAY_MUTATIONS_WHOLE_GENOME, "true");
        defaultValues.put(SHOW_SINGLE_TRACK_PANE_KEY, "false");
        defaultValues.put(PORT_ENABLED, "true");
        defaultValues.put(EXPAND_FEAUTRE_TRACKS, "false");
        defaultValues.put(SHOW_ATTRIBUTE_VIEWS_KEY, "true");
        defaultValues.put(SHOW_MISSING_DATA_KEY, "false");
        defaultValues.put(SHOW_SINGLE_TRACK_PANE_KEY, "false");
        defaultValues.put(SHOW_EXPAND_ICON, "false");
        defaultValues.put(SHOW_DEFAULT_TRACK_ATTRIBUTES, "false");

        defaultValues.put(IONTORRENT_FLOWDIST_HIDE_FIRST_HP, "true");
        defaultValues.put(IONTORRENT_FLOWDIST_BINSIZE, "15");
        defaultValues.put(IONTORRENT_FLOWDIST_CHARTTYPE, "LINE");
        defaultValues.put(IONTORRENT_SERVER, "ioneast.ite");
        defaultValues.put(IONTORRENT_RESULTS, "/results/analysis/output/Home/");

        defaultValues.put(CHART_DRAW_TOP_BORDER, "false");
        defaultValues.put(CHART_DRAW_BOTTOM_BORDER, "false");
        defaultValues.put(CHART_COLOR_BORDERS, "true");
        defaultValues.put(CHART_DRAW_TRACK_NAME, "false");
        defaultValues.put(CHART_DRAW_Y_AXIS, "false");
        defaultValues.put(CHART_AUTOSCALE, "false");
        defaultValues.put(CHART_SHOW_DATA_RANGE, "true");
        defaultValues.put(CHART_COLOR_TRACK_NAME, "true");
        defaultValues.put(CHART_TRACK_HEIGHT_KEY, "40");
        defaultValues.put(CHART_SHOW_ALL_HEATMAP, "false");
        defaultValues.put(UNLOAD_ON_GENOME_CHANGE, "false");

        defaultValues.put(SAM_SHOW_DUPLICATES, "false");
        defaultValues.put(SAM_SHOW_SOFT_CLIPPED, "false");
        defaultValues.put(SAM_FLAG_UNMAPPED_PAIR, "false");
        defaultValues.put(SAM_AUTO_SORT, "false");
        defaultValues.put(SAM_SHADE_CENTER, "true");
        defaultValues.put(SAM_SHOW_REF_SEQ, "false");
        defaultValues.put(SAM_SHOW_CENTER_LINE, "true");
        defaultValues.put(SAM_SHOW_COV_TRACK, "true");
        defaultValues.put(SAM_SHADE_BASES, ShadeBasesOption.QUALITY.toString());
        defaultValues.put(SAM_FILTER_ALIGNMENTS, "false");
        defaultValues.put(SAM_FILTER_SECONDARY_ALIGNMENTS, "false");
        defaultValues.put(SAM_FILTER_SUPPLEMENTARY_ALIGNMENTS, "false");
        defaultValues.put(SAM_FILTER_FAILED_READS, "true");
        defaultValues.put(SAM_DOWNSAMPLE_READS, "true");
        defaultValues.put(SAM_SAMPLING_WINDOW, "50");
        defaultValues.put(SAM_SAMPLING_COUNT, "100");
        defaultValues.put(SAM_BASE_QUALITY_MIN, "5");
        defaultValues.put(SAM_BASE_QUALITY_MAX, "20");
        defaultValues.put(SAM_FILTER_URL, null);
        defaultValues.put(SAM_QUALITY_THRESHOLD, "0");
        defaultValues.put(SAM_ALLELE_THRESHOLD, "0.2f");
        defaultValues.put(SAM_ALLELE_USE_QUALITY, "true");
        defaultValues.put(SAM_MIN_INSERT_SIZE_THRESHOLD, "50");
        defaultValues.put(SAM_MAX_INSERT_SIZE_THRESHOLD, "1000");
        defaultValues.put(SAM_MIN_INSERT_SIZE_PERCENTILE, "0.5");
        defaultValues.put(SAM_MAX_INSERT_SIZE_PERCENTILE, "99.5");
        defaultValues.put(SAM_MAX_VISIBLE_RANGE, "30");
        defaultValues.put(SAM_COLOR_BY, "UNEXPECTED_PAIR");
        defaultValues.put(SAM_COLOR_BY_TAG, "");
        defaultValues.put(SAM_GROUP_BY_TAG, "");
        defaultValues.put(SAM_SORT_BY_TAG, "");
        defaultValues.put(SAM_BISULFITE_CONTEXT, "CG");
        defaultValues.put(SAM_COMPUTE_ISIZES, "true");
        defaultValues.put(SAM_FLAG_ZERO_QUALITY, "true");
        defaultValues.put(SAM_SHOW_JUNCTION_TRACK, "false");
        defaultValues.put(SAM_JUNCTION_MIN_FLANKING_WIDTH, "0");
        defaultValues.put(SAM_JUNCTION_MIN_COVERAGE, "1");
        defaultValues.put(SAM_SHOW_JUNCTION_FLANKINGREGIONS, "true");
        defaultValues.put(SAM_NOMESEQ_ENABLED, "false");
        defaultValues.put(SAM_COUNT_DELETED_BASES_COVERED, "false");
        defaultValues.put(SAM_FLAG_LARGE_INSERTIONS, "false");
        defaultValues.put(SAM_LARGE_INSERTIONS_THRESHOLD, "1");
        defaultValues.put(SAM_SORT_OPTION, "NUCLEOTIDE");
        defaultValues.put(SAM_SHOW_GROUP_SEPARATOR, "true");
        defaultValues.put(SAM_COMPLETE_READS_ONLY, "false");
        defaultValues.put(SAM_SHOW_ALL_BASES, "false");

        defaultValues.put(SAM_REDUCED_MEMORY_MODE, "false");
        defaultValues.put(SAM_MIN_INDEL_SIZE, "25");

        defaultValues.put(SAM_SHOW_ALIGNMENT_TRACK, "true");

        defaultValues.put(BYPASS_FILE_AUTO_DISCOVERY, "false");

        defaultValues.put(NORMALIZE_COVERAGE, "false");

        defaultValues.put(SHOW_GENOME_SERVER_WARNING, "true");

        defaultValues.put(SEARCH_ZOOM, "true");


        defaultValues.put(GENOMES_SERVER_URL, DEFAULT_GENOME_URL);
        defaultValues.put(OVERLAY_ATTRIBUTE_KEY, "LINKING_ID");
        defaultValues.put(DEFAULT_GENOME_KEY, DEFAULT_GENOME);

        defaultValues.put(USE_PROXY, "false");
        defaultValues.put(PROXY_AUTHENTICATE, "false");
        defaultValues.put(PORT_NUMBER, "60151");
        defaultValues.put(TRACK_HEIGHT_KEY, "15");
        defaultValues.put(FLANKING_REGION, "2000");

        defaultValues.put(SHOW_SEQUENCE_TRANSLATION, "false");
        defaultValues.put(MAX_SEQUENCE_RESOLUTION, "2");

        defaultValues.put(AUTO_UPDATE_GENOMES, "true");

        defaultValues.put(GWAS_TRACK_HEIGHT, "200");
        defaultValues.put(GWAS_DESCRIPTION_CACHE_SIZE, "10000");
        defaultValues.put(GWAS_MIN_POINT_SIZE, "3");
        defaultValues.put(GWAS_MAX_POINT_SIZE, "7");
        defaultValues.put(GWAS_USE_CHR_COLORS, "true");
        defaultValues.put(GWAS_SINGLE_COLOR, "false");
        defaultValues.put(GWAS_ALTERNATING_COLORS, "false");
        defaultValues.put(GWAS_PRIMARY_COLOR, "69,101,183");
        defaultValues.put(GWAS_SECONDARY_COLOR, "250,169,10");
        defaultValues.put(GWAS_SHOW_AXIS, "true");

        defaultValues.put(DEFAULT_FONT_SIZE, "10");
        defaultValues.put(DEFAULT_FONT_FAMILY, "Arial");
        defaultValues.put(DEFAULT_FONT_ATTRIBUTE, String.valueOf(Font.PLAIN));
        defaultValues.put(SCALE_FONTS, "false");

        boolean isMac = System.getProperty("os.name").toLowerCase().startsWith("mac");
        defaultValues.put(ENABLE_ANTIALISING, String.valueOf(isMac));

        defaultValues.put(NAME_PANEL_WIDTH, "160");
        defaultValues.put(BACKGROUND_COLOR, "250,250,250");

        defaultValues.put(GENOME_SPACE_ENABLE, "true");
        defaultValues.put(GENOME_SPACE_DM_SERVER, "https://dm.genomespace.org/datamanager/v1.0/");
        defaultValues.put(GENOME_SPACE_ATM_SERVER, "https://atm.genomespace.org/atm/v1.0/");
        defaultValues.put(GENOME_SPACE_IDENTITY_SERVER, "https://identitydev.genomespace.org:8444/identityServer/basic");

        // Affective computing mode
        defaultValues.put(ENABLE_EXOME_BUTTON, "false");

        defaultValues.put(DB_ENABLED, "false");
        defaultValues.put(DB_HOST, "");
        defaultValues.put(DB_NAME, "");
        defaultValues.put(DB_PORT, "-1");

        String defaultDataURL = DEFAULT_DATA_URL;
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

        defaultValues.put(DATA_SERVER_URL_KEY, defaultDataURL);

        defaultValues.put(FRAME_STATE_KEY, "" + Frame.NORMAL);

        defaultValues.put(CBIO_MUTATION_THRESHOLD, "1");
        defaultValues.put(CBIO_AMPLIFICATION_THRESHOLD, "0.9");
        defaultValues.put(CBIO_DELETION_THRESHOLD, "0.9");
        defaultValues.put(CBIO_EXPRESSION_UP_THRESHOLD, "1.0");
        defaultValues.put(CBIO_EXPRESSION_DOWN_THRESHOLD, "1.0");

        defaultValues.put(TOOLTIP_INITIAL_DELAY, "50");
        defaultValues.put(TOOLTIP_RESHOW_DELAY, "50");
        defaultValues.put(TOOLTIP_DISMISS_DELAY, "60000");
        defaultValues.put(DETAILS_BEHAVIOR_KEY, IGVCommandBar.SHOW_DETAILS_BEHAVIOR.HOVER.name());

        defaultValues.put(SHOW_SIZE_WARNING, "true");

        defaultValues.put(SKIP_VERSION, "");

        defaultValues.put(COLOR_A, "0,150,0");
        defaultValues.put(COLOR_C, "0,0,255");
        defaultValues.put(COLOR_T, "255,0,0");
        defaultValues.put(COLOR_G, "209,113,5");
        defaultValues.put(COLOR_N, ColorUtilities.colorToString(Color.gray));
        defaultValues.put(SAM_COLOR_A, "0,255,0");
        defaultValues.put(SAM_COLOR_C, "0,0,255");
        defaultValues.put(SAM_COLOR_T, "255,0,0");
        defaultValues.put(SAM_COLOR_G, "209,113,5");
        defaultValues.put(SAM_COLOR_N, ColorUtilities.colorToString(Color.gray.brighter()));

        defaultValues.put(HOMREF_COLOR, "235,235,235");
        defaultValues.put(HETVAR_COLOR, "0,0,255");
        defaultValues.put(HOMVAR_COLOR, "0,245,255");
        defaultValues.put(NOCALL_COLOR, "255,255,255");
        defaultValues.put(AF_REF_COLOR, "0,0,220");
        defaultValues.put(AF_VAR_COLOR, "255,0,0");

        defaultValues.put(VARIANT_COLOR_BY_ALLELE_FREQ, "true");

        defaultValues.put(SASHIMI_SHOW_COVERAGE, "true");

        defaultValues.put(ENABLE_GOOGLE_MENU, "false");
        defaultValues.put(SAVE_GOOGLE_CREDENTIALS, "true");

        defaultValues.put(DEFAULT_VISIBILITY_WINDOW, "-1");

        defaultValues.put(BLAT_URL, "http://genome.cse.ucsc.edu/cgi-bin/hgBlat");

        defaultValues.put(GENE_LIST_BED_FORMAT, "false");

        defaultValues.put(SESSION_RELATIVE_PATH, "false");

        defaultValues.put(SHOW_LOS, "true");

    }

    /**
     * Set an alternate preference file. Mainly for testing.
     *
     * @param s
     */
    public void setPrefsFile(String s) {
        preferences = new IGVPreferences(new File(s));
        clearCaches();
    }

    public static String generateGenomeIdString(Collection<GenomeListItem> genomeListItems) {
        String genomeString = "";

        for (GenomeListItem serverItem : genomeListItems) {
            genomeString += serverItem.getId() + HISTORY_DELIMITER;
        }

        genomeString = genomeString.substring(0, genomeString.length() - 1);
        return genomeString;
    }

    public void saveGenomeIdDisplayList(Collection<GenomeListItem> genomeListItems) {
        preferences.put(GENOME_ID_DISPLAY_LIST_KEY, generateGenomeIdString(genomeListItems));
    }

    /**
     * Returns a list of genomeIds to display in the combo box
     *
     * @return
     */
    public String[] getGenomeIdDisplayList() {
        return getArray(GENOME_ID_DISPLAY_LIST_KEY);
    }

    /**
     * Get a property which is a delimited list of entries
     *
     * @param key
     * @return The string array of tokens, or an empty array if not present
     */
    private String[] getArray(String key) {
        String stringProp = get(key);
        if (stringProp == null) {
            return new String[0];
        } else {
            return stringProp.split(HISTORY_DELIMITER);
        }
    }

    /**
     * Get the path to the CLI plugin specified by the
     * Id and tool name.
     *
     * @param pluginId
     * @param toolName
     * @return
     * @see #putToolPath(String, String, String)
     * @see #genToolKey
     */
    public String getToolPath(String pluginId, String toolName) {
        return get(genToolKey(pluginId, toolName, "path"));
    }

    /**
     * Set the path to the CLI plugin
     *
     * @param pluginId
     * @param toolName
     * @param path
     * @see #getToolPath(String, String)
     * @see #genToolKey
     */
    public void putToolPath(String pluginId, String toolName, String path) {
        put(genToolKey(pluginId, toolName, "path"), path);
    }

    /**
     * Used to generate a unique string based on pluginId, toolName, and attribute
     *
     * @param pluginId
     * @param toolName
     * @param key
     * @return
     */
    private String genToolKey(String pluginId, String toolName, String key) {
        return String.format("%s:%s:%s", pluginId, toolName.replace(' ', '_'), key.replace(' ', '_'));
    }

    public void putArgumentValue(String pluginId, String toolName, String command, String argName, String argValue) {
        String key = genArgKey(pluginId, toolName, command, argName);
        put(key, argValue);
    }

    public String getArgumentValue(String pluginId, String toolName, String commandName, String argId) {
        return get(genArgKey(pluginId, toolName, commandName, argId));
    }

    private String genArgKey(String pluginId, String toolName, String command, String argId) {
        return genToolKey(pluginId, toolName, String.format("%s:%s", command, argId));
    }

    public String[] getIGVPluginList() {
        return getArray(IGV_PLUGIN_LIST_KEY);
    }

    /**
     * Returns a preference value from either the preferences,
     * or system property. If the value was provided as a system property,
     * is is saved to the preferences.
     * Intended usecase is features only usable by certain groups but
     * not intended for all of IGV.
     * <p/>
     * Example
     * java -jar Denable.tools=true igv.jar
     * <p/>
     * getPersistent("enable.tools", null) returns "true"
     * and saves it to preferences, so a subsequent invocation
     * will return true as well:
     * java -jar igv.jar
     * getPersistent("enable.tools", "false") returns "true" also
     *
     * @param key
     * @param def default value. NOT SAVED
     * @return
     * @see org.broad.igv.session.Session#getPersistent(String, String)
     */
    public String getPersistent(String key, String def) {
        String value = System.getProperty(key);
        if (value != null) {
            put(key, value);
            return value;
        } else {
            return get(key, def);
        }
    }


    /**
     * List of keys that affect the alignments loaded.  This list is used to trigger a reload, if required.
     * Not all alignment preferences need trigger a reload, this is a subset.
     */
    static java.util.List<String> SAM_RELOAD_KEYS = Arrays.asList(
            PreferenceManager.SAM_QUALITY_THRESHOLD,
            PreferenceManager.SAM_FILTER_ALIGNMENTS,
            PreferenceManager.SAM_FILTER_URL,
            PreferenceManager.SAM_MAX_VISIBLE_RANGE,
            PreferenceManager.SAM_SHOW_DUPLICATES,
            PreferenceManager.SAM_SHOW_SOFT_CLIPPED,
            PreferenceManager.SAM_SAMPLING_COUNT,
            PreferenceManager.SAM_SAMPLING_WINDOW,
            PreferenceManager.SAM_FILTER_FAILED_READS,
            PreferenceManager.SAM_DOWNSAMPLE_READS,
            PreferenceManager.SAM_FILTER_SECONDARY_ALIGNMENTS,
            PreferenceManager.SAM_FILTER_SUPPLEMENTARY_ALIGNMENTS,
            PreferenceManager.SAM_JUNCTION_MIN_FLANKING_WIDTH,
            PreferenceManager.SAM_JUNCTION_MIN_COVERAGE
    );


}
