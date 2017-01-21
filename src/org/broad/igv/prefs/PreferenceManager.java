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
package org.broad.igv.prefs;


import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.feature.genome.GenomeListItem;
import org.broad.igv.renderer.ColorScaleFactory;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.sam.AlignmentTrack.ShadeBasesOption;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.AboutDialog;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.IGVCommandBar;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.color.PaletteColorTable;
import org.broad.igv.ui.event.AlignmentTrackEvent;
import org.broad.igv.ui.util.PropertyManager;
import org.broad.igv.util.HttpUtils;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.util.*;

import static org.broad.igv.prefs.Constants.*;

/**
 * Manages user preferences.
 */
public class PreferenceManager implements PropertyManager {

    private static Logger log = Logger.getLogger(PreferenceManager.class);

    IGVPreferences preferences;

    Map<String, String> defaultValues;


    // Cached preference values
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

            boolean refreshSAM = false;
            for (String key: SAM_REFRESH_KEYS) {
                if (updatedPreferenceMap.containsKey(key)) {
                    refreshSAM = true;
                    break;
                }
            }

            if (reloadSAM) {
                if (updatedPreferenceMap.containsKey(SAM_MAX_VISIBLE_RANGE)) {
                    igv.notifyAlignmentTrackEvent(this, AlignmentTrackEvent.Type.VISIBILITY_WINDOW);
                }
                igv.notifyAlignmentTrackEvent(this, AlignmentTrackEvent.Type.RELOAD);
            }
            // A reload is harsher than a refresh; only send the weaker request if the stronger one is not sent.
            if (!reloadSAM && refreshSAM) {
                igv.notifyAlignmentTrackEvent(this, AlignmentTrackEvent.Type.REFRESH);
            }
            if (updatedPreferenceMap.containsKey(SAM_ALLELE_THRESHOLD)) {
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
        return get(GENOMES_SERVER_URL);
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
        put(SHOW_ATTRIBUTE_VIEWS_KEY, Boolean.toString(isShowable));
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
    public static ContinuousColorScale getDefaultColorScale(TrackType type) {
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
                cs = getDefaultColorScale(Color.BLUE, Color.WHITE, Color.RED);
                cs.setNoDataColor(new Color(225, 225, 225));
                return cs;

            case COPY_NUMBER:
            case ALLELE_SPECIFIC_COPY_NUMBER:
            case CNV:
                return getDefaultColorScale(Color.BLUE, Color.WHITE, Color.RED);

            default:
                return null;
        }
    }

    public static ContinuousColorScale getDefaultColorScale(Color negColor, Color neutralColor, Color posColor) {
        return new ContinuousColorScale(-0.1, -1.5, 0.1, 1.5, negColor, neutralColor, posColor);
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

        remove(USE_PROXY);
        remove(PROXY_HOST);
        remove(PROXY_PORT);
        remove(PROXY_AUTHENTICATE);
        remove(PROXY_USER);
        remove(PROXY_PW);
        remove(PROXY_TYPE);
        remove(PROXY_WHITELIST);
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
        defaultValues.put(SAM_QUICK_CONSENSUS_MODE, "false");
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
        defaultValues.put(SAM_HIDDEN_TAGS, "SA,MD,XA,RG");
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
        defaultValues.put(SAM_GROUP_BY_POS, "");
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
        defaultValues.put(SAM_FLAG_LARGE_INDELS, "true");
        defaultValues.put(SAM_LARGE_INDELS_THRESHOLD, "1");
        defaultValues.put(SAM_FLAG_CLIPPING, "false");
        defaultValues.put(SAM_CLIPPING_THRESHOLD, "0");
        defaultValues.put(SAM_SORT_OPTION, "NUCLEOTIDE");
        defaultValues.put(SAM_GROUP_OPTION, "NONE");
        defaultValues.put(SAM_SHOW_GROUP_SEPARATOR, "true");
        defaultValues.put(SAM_COMPLETE_READS_ONLY, "false");
        defaultValues.put(SAM_SHOW_ALL_BASES, "false");

        defaultValues.put(SAM_REDUCED_MEMORY_MODE, "false");

        defaultValues.put(SAM_HIDE_SMALL_INDEL_BP, "false");
        defaultValues.put(SAM_SMALL_INDEL_BP_THRESHOLD, "0");

        defaultValues.put(SAM_LINK_READS, "false");
        defaultValues.put(SAM_LINK_TAG, "READNAME");

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
            SAM_QUALITY_THRESHOLD,
            SAM_FILTER_ALIGNMENTS,
            SAM_FILTER_URL,
            SAM_MAX_VISIBLE_RANGE,
            SAM_SHOW_DUPLICATES,
            SAM_SHOW_SOFT_CLIPPED,
            SAM_SAMPLING_COUNT,
            SAM_SAMPLING_WINDOW,
            SAM_FILTER_FAILED_READS,
            SAM_DOWNSAMPLE_READS,
            SAM_FILTER_SECONDARY_ALIGNMENTS,
            SAM_FILTER_SUPPLEMENTARY_ALIGNMENTS,
            SAM_JUNCTION_MIN_FLANKING_WIDTH,
            SAM_JUNCTION_MIN_COVERAGE
    );

    /**
     * List of keys that do not affect the alignments loaded but do affect how those
     * alignments are drawn.  A refresh is softer than a reload.
    */
    static java.util.List<String> SAM_REFRESH_KEYS = Arrays.asList(
        SAM_QUICK_CONSENSUS_MODE,
        SAM_ALLELE_THRESHOLD,
        SAM_FLAG_LARGE_INDELS,
        SAM_LARGE_INDELS_THRESHOLD
    );

}
