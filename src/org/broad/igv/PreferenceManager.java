/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.maf.MAFManager;
import org.broad.igv.renderer.ColorScaleFactory;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.AboutDialog;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.util.PropertyManager;
import org.broad.igv.ui.util.ColorTable;
import org.broad.igv.util.IGVHttpUtils;

import static org.broad.igv.ui.util.UIUtilities.getcommaSeparatedRGBString;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.List;

/**
 * @author eflakes
 */
public class PreferenceManager implements PropertyManager {

    private static Logger log = Logger.getLogger(PreferenceManager.class);

    private static final String DEFAULT_PORT_NUMBER = "60151";
    private static final String DEFAULT_AUTOSCALE = "true";
    private static final String DEFAULT_SHOW_REFERENCE = "false";

    public static final String INITIAL_TRACK_HEIGHT = "15";
    private final String SYSTEM_DEFAULT_FOR_DIRECT_DRAW = System.getProperty("sun.java2d.noddraw");
    /**
     * Property keys below
     */
    public static final String SAM_SHOW_ALL_MISMATCHES = "COVERAGE.SHOW_ALL_MISMATCHES";
    public static final String COVERAGE_SHOW_ALL_MISMATCHES = SAM_SHOW_ALL_MISMATCHES;
    public static final String CHART_DRAW_TOP_BORDER = "CHART.DRAW_TOP_BORDER";
    public static final String CHART_DRAW_BOTTOM_BORDER = "CHART.DRAW_BOTTOM_BORDER";
    public static final String CHART_COLOR_BORDERS = "CHART.COLOR_BORDERS";
    public static final String CHART_DRAW_Y_AXIS = "CHART.DRAW_AXIS";
    public static final String CHART_DRAW_TRACK_NAME = "CHART.DRAW_TRACK_NAME";
    public static final String CHART_COLOR_TRACK_NAME = "CHART.COLOR_TRACK_NAME";
    public static final String CHART_AUTOSCALE = "CHART.AUTOSCALE";
    public static final String CHART_SHOW_DATA_RANGE = "CHART.SHOW_DATA_RANGE";

    public static final String SAM_ALLELE_THRESHOLD = "SAM.ALLELE_THRESHOLD";
    public static final String SAM_QUALITY_THRESHOLD = "SAM.QUALITY_THRESHOLD";
    public static final String SAM_MAX_INSERT_SIZE_THRESHOLD = "SAM.INSERT_SIZE_THRESHOLD";
    public static final String SAM_MIN_INSERT_SIZE_THRESHOLD = "SAM.MIN_INSERT_SIZE_THRESHOLD";
    public static final String SAM_MAX_INSERT_SIZE_PERCENTILE = "SAM.ISIZE_MAX_PERCENTILE";
    public static final String SAM_MIN_INSERT_SIZE_PERCENTILE= "SAM.MIN_ISIZE_MIN_PERCENTILE";
    public static final String SAM_COMPUTE_INSERT_SIZE_THRESHOLD = "SAM.COMPUTE_ISZIE";
    public static final String SAM_AUTO_SORT = "SAM.AUTOSORT";
    public static final String SAM_SHADE_CENTER = "SAM.SHADE_CENTER";
    public static final String SAM_SHOW_REF_SEQ = "SAM.SHOW_REF_SEQ";
    public static final String SAM_SHOW_COV_TRACK = "SAM.SHOW_COV_TRACK";
    public static final String SAM_SHADE_BASE_QUALITY = "SAM.SHADE_BASE_QUALITY";
    public static final String SAM_BASE_QUALITY_MIN = "SAM.BASE_QUALITY_MIN";
    public static final String SAM_BASE_QUALITY_MAX = "SAM.BASE_QUALITY_MAX";
    public static final String SAM_FILTER_ALIGNMENTS = "SAM.FILTER_ALIGNMENTS";
    public static final String SAM_FILTER_URL = "SAM.FILTER_URL";
    public static final String SAM_MAX_VISIBLE_RANGE = "SAM.MAX_VISIBLE_RANGE";
    public static final String SAM_SHOW_DUPLICATES = "SAM.SHOW_DUPLICATES";
    public static final String SAM_SHOW_SOFT_CLIPPED = "SAM.SHOW_SOFT_CLIPPED";
    public static final String SAM_FLAG_UNMAPPED_PAIR = "SAM.FLAG_UNMAPPED_PAIR";
    public static final String SAM_MAX_LEVELS = "SAM.MAX_LEVELS";
    public static final String SAM_MAX_READS = "SAM.MAX_READS";
    public static final String SAM_COLOR_BY = "SAM.COLOR_BY";
    public static final String SAM_FILTER_FAILED_READS = "SAM.FILTER_FAILED_READS";
    public static final String SAM_COMPUTE_ISIZES = "SAM.COMPUTE_ISIZES";
    //dhmay adding 20110208
    public static final String SAM_SHOW_JUNCTION_TRACK = "SAM.SHOW_JUNCTION_TRACK";


    public static final String EXPAND_FEAUTRE_TRACKS = "EXPAND_FEATURE_TRACKS";
    public static final String PORT_ENABLED = "PORT_ENABLED";
    public static final String PORT_NUMBER = "PORT_NUMBER";
    public static final String COLOR_SCALE_KEY = "COLOR_SCALE_";
    final public static String FRAME_BOUNDS_KEY = "IGV.Bounds";
    final public static String RECENT_SESSION_KEY = "IGV.Session.recent.sessions";
    final public static String TRACK_HEIGHT_KEY = "IGV.track.height";
    final public static String CHART_TRACK_HEIGHT_KEY = "IGV.chart.track.height";
    final public static String SHOW_MISSING_DATA_KEY = "IGV.track.show.missing.data";
    final public static String SHOW_ATTRIBUTE_VIEWS_KEY = "IGV.track.show.attribute.views";
    final public static String SHOW_SINGLE_TRACK_PANE_KEY = "IGV.single.track.pane";
    final public static String GENOMES_SERVER_URL = "IGV.genome.sequence.dir";
    final public static String JOIN_ADJACENT_SEGMENTS_KEY = "IGV.join.adjacent.segments";
    final public static String SHOW_REGION_BARS = "SHOW_REGION_BARS";
    final public static String LAST_EXPORTED_REGION_DIRECTORY = "LAST_EXPORTED_REGION_DIRECTORY";
    final static public String LAST_TRACK_DIRECTORY = "LAST_TRACK_DIRECTORY";
    final static public String LAST_SNAPSHOT_DIRECTORY = "LAST_SNAPSHOT_DIRECTORY";
    final static public String LAST_GENOME_IMPORT_DIRECTORY = "LAST_GENOME_IMPORT_DIRECTORY";
    final static public String LAST_SESSION_DIRECTORY = "LAST_SESSION_DIRECTORY";
    final static public String DEFAULT_GENOME_KEY = "DEFAULT_GENOME_KEY";
    final static public String LAST_CHROMOSOME_VIEWED_KEY = "LAST_CHROMOSOME_VIEWED_KEY";
    final public static String MUTATION_INDEL_COLOR_KEY = "MUTATION_INDEL_COLOR_KEY";
    final public static String MUTATION_MISSENSE_COLOR_KEY = "MUTATION_MISSENSE_COLOR_KEY";
    final public static String MUTATION_NONSENSE_COLOR_KEY = "MUTATION_NONSENSE_COLOR_KEY";
    final public static String MUTATION_SPLICE_SITE_COLOR_KEY = "MUTATION_SPLICE_SITE_COLOR_KEY";
    final public static String MUTATION_SYNONYMOUS_COLOR_KEY = "MUTATION_SYNONYMOUS_COLOR_KEY";
    final public static String MUTATION_TARGETED_REGION_COLOR_KEY = "MUTATION_TARGETED_REGION_COLOR_KEY";
    final public static String MUTATION_UNKNOWN_COLOR_KEY = "MUTATION_UNKNOWN_COLOR_KEY";
    final public static String OVERLAY_TRACKS_KEY = "OVERLAY_TRACKS_KEY";
    final public static String DISPLAY_OVERLAY_TRACKS_KEY = "DISPLAY_OVERLAY_TRACKS_KEY";
    final public static String OVERLAY_ATTRIBUTE_KEY = "OVERLAY_ATTRIBUTE_KEY";
    final public static String COLOR_OVERLAY_KEY = "COVER_OVERLAY_KEY";
    final public static String ENABLE_LINKED_SORTING = "ENABLE_LINKED_SORTING";
    final public static String TRACK_ATTRIBUTE_NAME_KEY = "TRACK_ATTRIBUTE_NAME_KEY";
    final public static String DATA_SERVER_URL_KEY = "MASTER_RESOURCE_FILE_KEY";
    final public static String CHECKED_RESOURCES_KEY = "CHECKED_RESOURCES_KEY";
    final public static String DEFINE_GENOME_INPUT_DIRECTORY_KEY = "DEFINE_GENOME_INPUT_DIRECTORY_KEY";
    final public static String LAST_CYTOBAND_DIRECTORY_KEY = "LAST_CYTOBAND_DIRECTORY_KEY";
    final public static String LAST_REFFLAT_DIRECTORY_KEY = "LAST_REFFLAT_DIRECTORY_KEY";
    final public static String LAST_FASTA_DIRECTORY_KEY = "LAST_FASTA_DIRECTORY_KEY";
    final public static String LAST_SEQUENCE_DIRECTORY_KEY = "LAST_SEQUENCE_DIRECTORY_KEY";
    final public static String MAF_SPECIES_KEY = "MAF_SPECIES_KEY";
    final public static String PROBE_MAPPING_KEY = "PROBE_MAPPING_KEY";
    final public static String SEARCH_ZOOM = "SEARCH_ZOOM";
    final public static String NORMALIZE_COVERAGE = "NORMALIZE_COVERAGE";
    public static final String SHOW_EXPAND_ICON = "SHOW_EXPAND_ICON";
    public static final String FTP_ANON_EMAIL = "FTP_ANON_EMAIL";

    public static final String SHOW_SIZE_WARNING = "SHOW_SIZE_WARNING";
    public static final String SHOW_GENOME_SERVER_WARNING = "SHOW_GENOME_SERVER_WARNING";

    final public static String USE_PROXY = "PROXY.USE";
    final public static String PROXY_HOST = "PROXY.HOST";
    final public static String PROXY_PORT = "PROXY.PORT";
    final public static String PROXY_AUTHENTICATE = "PROXY.AUTHENTICATE";
    final public static String PROXY_USER = "PROXY.USERNAME";
    final public static String PROXY_PW = "PROXY.PW";

    final public static String KNOWN_SNPS = "KNOWN_SNPS_FILE";

    final public static String USE_BYTE_RANGE = "UseHttpByteRange";

    public static final String FLANKING_REGION = "FLAKING_REGIONS";

    public static final String SHOW_SEQUENCE_TRANSLATION = "SHOW_SEQUENCE_TRANSLATION";
    public static final String MAX_SEQUENCE_RESOLUTION = "MAX_SEQUENCE_RESOLUTION";

    public static final String AUTO_UPDATE_GENOMES = "AUTO_UPDATE_GENOMES";

    final public static String GWAS_TRACK_HEIGHT = "GWAS_TRACK_HEIGHT";
    final public static String GWAS_DESCRIPTION_CACHE_SIZE = "GWAS_DESCRIPTION_CACHE_SIZE";
    final public static String GWAS_MIN_POINT_SIZE = "GWAS_MIN_POINT_SIZE";
    final public static String GWAS_MAX_POINT_SIZE = "GWAS_MAX_POINT_SIZE";
    final public static String GWAS_USE_CHR_COLORS = "GWAS_USE_CHR_COLORS";
    final public static String GWAS_SINGLE_COLOR = "GWAS_SINGLE_COLOR";
    final public static String GWAS_ALTERNATING_COLORS = "GWAS_ALTERNATING_COLORS";
    final public static String GWAS_PRIMARY_COLOR = "GWAS_PRIMARY_COLOR";
    final public static String GWAS_SECONDARY_COLOR = "GWAS_SECONDARY_COLOR";

    public static final String DEFAULT_FONT_SIZE = "DEFAULT_FONT_SIZE";
    public static final String DEFAULT_FONT_FAMILY = "DEFAULT_FONT_FAMILY";
    public static final String DEFAULT_FONT_ATTRIBUTE = "DEFAULT_FONT_ATTRIBUTE";

    public static final String NAME_PANEL_WIDTH = "NAME_PANEL_WIDTH";





    public static String DEFAULT_DATA_SERVER_URL;

    static {
        Properties properties = new Properties();
        try {
            properties.load(AboutDialog.class.getResourceAsStream("/resources/about.properties"));
            DEFAULT_DATA_SERVER_URL = properties.getProperty("master-resource-url", "http://www.broadinstitute.org/igvdata/$$_dataServerRegistry.txt");
            if (DEFAULT_DATA_SERVER_URL.equals("@DEFAULT_MASTER_RESOURCE_URL")) {
                DEFAULT_DATA_SERVER_URL = "http://www.broadinstitute.org/igvdata/$$_dataServerRegistry.txt";
            }
        } catch (IOException e) {
            DEFAULT_DATA_SERVER_URL = "http://www.broadinstitute.org/igvdata/$$_dataServerRegistry.txt";
        }
    }

    /**
     * The preference cache
     */
    IGVPreferences preferences;

    /**
     * Default values
     */
    Map<String, String> defaultValues;


    /**
     * Cache of preference values.  Profiling reveals that Preferences.get()
     * is taking huge amounts of time.  There are hundereds of thousands of
     * calls to this to get the track height,  this is possibly a bad design
     * decision, however caching the preference values solves the peformance
     * problem for now.
     */
    private Map<String, Boolean> booleanCache = new Hashtable();
    private Map<String, Object> objectCache = new Hashtable();
    private Map<TrackType, ContinuousColorScale> colorScaleCache = new Hashtable();

    public static PreferenceManager getInstance() {

        return instance;
    }

    private static PreferenceManager instance = new PreferenceManager();

    private PreferenceManager() {
        preferences = new IGVPreferences();
        initDefaultValues();
    }


    public String get(String key, String defaultString) {
        return preferences.get(key, defaultString);
    }

    public String get(String key) {
        return get(key, defaultValues.get(key));
    }


    /**
     * Return the preference as a boolean value.
     *
     * @param key
     * @return
     */
    public boolean getAsBoolean(String key) {
        Boolean boolValue = booleanCache.get(key);
        if (boolValue == null) {
            String value = get(key);
            boolValue = new Boolean(get(key, value));
            if (boolValue == null) {
                log.error("No default value for: " + key);
                return false;
            }
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
        Number value = (Number) objectCache.get(key);
        if (value == null) {
            String defValue = get(key);
            value = new Integer(get(key, defValue));
            if (value == null) {
                log.error("No default value for: " + key);
                return 0;
            }
            objectCache.put(key, value);
        }
        return value.intValue();
    }

    /**
     * Return the preference as an float.
     *
     * @param key
     * @return
     */
    public float getAsFloat(String key) {
        Number value = (Number) objectCache.get(key);
        if (value == null) {
            String defValue = get(key);
            value = new Float(get(key, defValue));
            if (value == null) {
                log.error("No default value for: " + key);
                return 0;
            }
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
        preferences.put(key, value);
        updateCaches(key, value);
    }

    public void put(String key, boolean b) {
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
        clearCaches();

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


    public List<String> getMafSpecies() {
        String tmp = get(MAF_SPECIES_KEY, null);

        String[] species = null;
        if (tmp == null) {
            species = MAFManager.species;
        } else {
            species = tmp.split(":");
        }
        return Arrays.asList(species);
    }

    public void setMafSpecies(List<String> species) {
        StringBuffer buf = new StringBuffer(species.size() * 7);
        Iterator<String> iter = species.iterator();
        while (iter.hasNext()) {
            buf.append(iter.next());
            if (iter.hasNext()) {
                buf.append(":");
            }
        }
        put(MAF_SPECIES_KEY, buf.toString());

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

        String lastFilePath = get(DEFINE_GENOME_INPUT_DIRECTORY_KEY, Globals.getUserDirectory().getAbsolutePath());

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

        String lastFilePath = get(LAST_GENOME_IMPORT_DIRECTORY, Globals.getUserDirectory().getAbsolutePath());

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
        String masterResourceFile = get(DATA_SERVER_URL_KEY, DEFAULT_DATA_SERVER_URL);
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

        String genome = get(DEFAULT_GENOME_KEY, Globals.DEFAULT_GENOME);
        return genome;
    }

    public void setLastTrackDirectory(File directory) {
        String lastDirectory = directory.getAbsolutePath();
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
                cs = new ContinuousColorScale(0.4, 1, Color.WHITE, Color.GREEN.darker());
                cs.setNoDataColor(new Color(225, 225, 225));
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
     * @param colorScheme
     */
    public void setMutationColorScheme(ColorTable colorScheme) {

        setMutationColorScheme(
                getcommaSeparatedRGBString(colorScheme.get("Indel")),
                getcommaSeparatedRGBString(colorScheme.get("Missense")),
                getcommaSeparatedRGBString(colorScheme.get("Nonsense")),
                getcommaSeparatedRGBString(colorScheme.get("Splice_site")),
                getcommaSeparatedRGBString(colorScheme.get("Synonymous")),
                getcommaSeparatedRGBString(colorScheme.get("Targeted_Region")),
                getcommaSeparatedRGBString(colorScheme.get("Unknown")));
    }

    /**
     * @param indel
     * @param missense
     * @param nonsense
     * @param spliceSite
     * @param synonymous
     * @param targetedRegion
     * @param unknown
     */
    public void setMutationColorScheme(String indel, String missense, String nonsense,
                                       String spliceSite, String synonymous, String targetedRegion,
                                       String unknown) {

        put(MUTATION_INDEL_COLOR_KEY, indel);
        put(MUTATION_MISSENSE_COLOR_KEY, missense);
        put(MUTATION_NONSENSE_COLOR_KEY, nonsense);
        put(MUTATION_SPLICE_SITE_COLOR_KEY, spliceSite);
        put(MUTATION_SYNONYMOUS_COLOR_KEY, synonymous);
        put(MUTATION_TARGETED_REGION_COLOR_KEY, targetedRegion);
        put(MUTATION_UNKNOWN_COLOR_KEY, unknown);
    }

    /**
     * @return
     */
    public ColorTable getMutationColorScheme() {

        String indelColor = get(MUTATION_INDEL_COLOR_KEY, null);
        String missenseColor = get(MUTATION_MISSENSE_COLOR_KEY, null);
        String nonsenseColor = get(MUTATION_NONSENSE_COLOR_KEY, null);
        String spliceSiteColor = get(MUTATION_SPLICE_SITE_COLOR_KEY, null);
        String synonymousColor = get(MUTATION_SYNONYMOUS_COLOR_KEY, null);
        String targetedRegionColor = get(MUTATION_TARGETED_REGION_COLOR_KEY, null);
        String unknownColor = get(MUTATION_UNKNOWN_COLOR_KEY, null);

        ColorTable colorTable = new ColorTable();
        if ((indelColor != null) && (missenseColor != null) && (nonsenseColor != null) &&
                (spliceSiteColor != null) &&
                (synonymousColor != null) &&
                (targetedRegionColor != null) &&
                (unknownColor != null)) {


            String rgb[] = indelColor.split(",");
            Color color1 = new Color(Integer.parseInt(rgb[0]), Integer.parseInt(rgb[1]), Integer.parseInt(rgb[2]));
            colorTable.put("Indel", color1);

            rgb = missenseColor.split(",");
            Color color2 = new Color(Integer.parseInt(rgb[0]), Integer.parseInt(rgb[1]), Integer.parseInt(rgb[2]));
            colorTable.put("Missense", color2);

            rgb = nonsenseColor.split(",");
            Color color3 = new Color(Integer.parseInt(rgb[0]), Integer.parseInt(rgb[1]), Integer.parseInt(rgb[2]));
            colorTable.put("Nonsense", color3);

            rgb = spliceSiteColor.split(",");
            Color color4 = new Color(Integer.parseInt(rgb[0]), Integer.parseInt(rgb[1]), Integer.parseInt(rgb[2]));
            colorTable.put("Splice_site", color4);

            rgb = synonymousColor.split(",");
            Color color5 = new Color(Integer.parseInt(rgb[0]), Integer.parseInt(rgb[1]),
                    Integer.parseInt(rgb[2]));
            colorTable.put("Synonymous", color5);

            rgb = targetedRegionColor.split(",");
            Color color6 = new Color(Integer.parseInt(rgb[0]), Integer.parseInt(rgb[1]),
                    Integer.parseInt(rgb[2]));
            colorTable.put("Targeted_Region", color6);

            rgb = unknownColor.split(",");
            Color color7 = new Color(Integer.parseInt(rgb[0]), Integer.parseInt(rgb[1]),
                    Integer.parseInt(rgb[2]));
            colorTable.put("Unknown", color7);
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
        IGVHttpUtils.updateProxySettings();
    }


    private void initDefaultValues() {

        defaultValues = new HashMap();

        defaultValues.put(PROBE_MAPPING_KEY, "false");
        defaultValues.put(SHOW_REGION_BARS, "false");
        defaultValues.put(JOIN_ADJACENT_SEGMENTS_KEY, "false");
        defaultValues.put(USE_BYTE_RANGE, "true");
        defaultValues.put(OVERLAY_TRACKS_KEY, "true");
        defaultValues.put(COLOR_OVERLAY_KEY, "false");
        defaultValues.put(ENABLE_LINKED_SORTING, "true");
        defaultValues.put(SHOW_SINGLE_TRACK_PANE_KEY, "false");
        defaultValues.put(PORT_ENABLED, "true");
        defaultValues.put(EXPAND_FEAUTRE_TRACKS, "false");
        defaultValues.put(SHOW_ATTRIBUTE_VIEWS_KEY, "true");
        defaultValues.put(SHOW_MISSING_DATA_KEY, "false");
        defaultValues.put(SHOW_SINGLE_TRACK_PANE_KEY, "false");
        defaultValues.put(SHOW_EXPAND_ICON, "false");
        defaultValues.put(DISPLAY_OVERLAY_TRACKS_KEY, "true");

        defaultValues.put(CHART_DRAW_TOP_BORDER, "false");
        defaultValues.put(CHART_DRAW_BOTTOM_BORDER, "false");
        defaultValues.put(CHART_COLOR_BORDERS, "true");
        defaultValues.put(CHART_DRAW_TRACK_NAME, "false");
        defaultValues.put(CHART_DRAW_Y_AXIS, "false");
        defaultValues.put(CHART_AUTOSCALE, "true");
        defaultValues.put(CHART_SHOW_DATA_RANGE, "true");

        defaultValues.put(CHART_COLOR_TRACK_NAME, "true");
        defaultValues.put(CHART_TRACK_HEIGHT_KEY, "40");

        defaultValues.put(SAM_SHOW_DUPLICATES, "false");
        defaultValues.put(SAM_SHOW_SOFT_CLIPPED, "false");
        defaultValues.put(SAM_FLAG_UNMAPPED_PAIR, "false");
        defaultValues.put(SAM_AUTO_SORT, "false");
        defaultValues.put(SAM_SHADE_CENTER, "true");
        defaultValues.put(SAM_SHOW_REF_SEQ, "false");
        defaultValues.put(SAM_SHOW_COV_TRACK, "true");
        defaultValues.put(SAM_SHADE_BASE_QUALITY, "true");
        defaultValues.put(SAM_FILTER_ALIGNMENTS, "false");
        defaultValues.put(SAM_FILTER_FAILED_READS, "true");
        defaultValues.put(SAM_MAX_LEVELS, "100");
        defaultValues.put(SAM_MAX_READS, "500000");
        defaultValues.put(SAM_BASE_QUALITY_MIN, "5");
        defaultValues.put(SAM_BASE_QUALITY_MAX, "20");
        defaultValues.put(SAM_SHOW_ALL_MISMATCHES, "false");
        defaultValues.put(SAM_FILTER_URL, null);
        defaultValues.put(SAM_QUALITY_THRESHOLD, "0");
        defaultValues.put(SAM_ALLELE_THRESHOLD, "0.2f");
        defaultValues.put(SAM_MIN_INSERT_SIZE_THRESHOLD, "50");
        defaultValues.put(SAM_MAX_INSERT_SIZE_THRESHOLD, "1000");
        defaultValues.put(SAM_MIN_INSERT_SIZE_PERCENTILE, "0.5");
        defaultValues.put(SAM_MAX_INSERT_SIZE_PERCENTILE, "99.5");
        defaultValues.put(SAM_COMPUTE_INSERT_SIZE_THRESHOLD, "false");
        defaultValues.put(SAM_MAX_VISIBLE_RANGE, "30");
        defaultValues.put(SAM_COLOR_BY, "INSERT_SIZE");
        defaultValues.put(SAM_COMPUTE_ISIZES, "false");
        defaultValues.put(SAM_SHOW_JUNCTION_TRACK, "false");


        defaultValues.put(NORMALIZE_COVERAGE, "false");

        defaultValues.put(SHOW_GENOME_SERVER_WARNING, "true");
        defaultValues.put(SHOW_SIZE_WARNING, "true");

        defaultValues.put(SEARCH_ZOOM, "true");


        defaultValues.put(PreferenceManager.GENOMES_SERVER_URL, UIConstants.DEFAULT_SERVER_GENOME_ARCHIVE_LIST);
        defaultValues.put(OVERLAY_ATTRIBUTE_KEY, "LINKING_ID");
        defaultValues.put(DEFAULT_GENOME_KEY, Globals.DEFAULT_GENOME);
        defaultValues.put(DATA_SERVER_URL_KEY, DEFAULT_DATA_SERVER_URL);

        defaultValues.put(USE_PROXY, "false");
        defaultValues.put(PROXY_AUTHENTICATE, "false");
        defaultValues.put(PORT_NUMBER, "60151");
        defaultValues.put(TRACK_HEIGHT_KEY, "15");
        defaultValues.put(OVERLAY_ATTRIBUTE_KEY, "LINKING_ID");

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

        defaultValues.put(DEFAULT_FONT_SIZE, "10");
        defaultValues.put(DEFAULT_FONT_FAMILY, "Arial");
        defaultValues.put(DEFAULT_FONT_ATTRIBUTE, String.valueOf(Font.PLAIN));

        defaultValues.put(NAME_PANEL_WIDTH, "160");

    }
}
