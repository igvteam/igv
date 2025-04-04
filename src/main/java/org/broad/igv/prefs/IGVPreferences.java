
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


import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.batch.CommandListener;
import org.broad.igv.event.AlignmentTrackEvent;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.feature.genome.GenomeListItem;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.renderer.ColorScaleFactory;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.renderer.SequenceRenderer;
import org.broad.igv.sam.mods.BaseModificationColors;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.*;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.color.PaletteColorTable;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.HttpUtils;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

import static org.broad.igv.prefs.Constants.*;

/**
 * Manages user preferences.
 */
public class IGVPreferences {

    private static Logger log = LogManager.getLogger(IGVPreferences.class);

    IGVPreferences parent;

    Map<String, String> userPreferences;
    Map<String, String> defaults;
    // Preferences which should persist for this session only
    Set<String> overrideKeys = new HashSet<>();

    // Cached non-string preference values
    private Map<String, Boolean> booleanCache = new Hashtable<>();
    private Map<String, Object> objectCache = new Hashtable<>();
    private Map<TrackType, ContinuousColorScale> colorScaleCache = new Hashtable<>();
    private PaletteColorTable mutationColorScheme = null;

    public IGVPreferences() {
        this.userPreferences = new HashMap<>();
    }

    public IGVPreferences(Map<String, String> userPreferences,
                          Map<String, String> defaults,
                          IGVPreferences parent) {
        this.parent = parent;
        this.defaults = defaults;
        this.userPreferences = userPreferences == null ? new HashMap<>() : userPreferences;
    }

    public String get(String key) {
        key = key.trim();

        if (userPreferences.containsKey(key)) {
            return userPreferences.get(key);
        } else if (parent != null && parent.userPreferences.containsKey(key)) {
            return parent.userPreferences.get(key);
        } else if (defaults != null && defaults.containsKey(key)) {
            return defaults.get(key);
        } else if (parent != null) {
            return parent.get(key);
        } else {
            return null;
        }
    }


    /**
     * Return preference with given key and specified default value.  If key is not present defaultValue is returned,
     * no search through defaults or hierarchy is performed.
     *
     * @param key
     * @param defaultValue
     * @return
     */
    public String get(String key, String defaultValue) {
        key = key.trim();
        String val = userPreferences.get(key);
        return val == null ? defaultValue : val;
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
                log.warn("No default value for: " + key);
                return false;
            }
            boolValue = Boolean.valueOf(get(key, value));
            booleanCache.put(key, boolValue);
        }
        return boolValue;
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
                log.warn("No default value for: " + key);
                return 0;
            }
            value = Integer.valueOf(get(key, defValue));
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
                log.warn("No default value for: " + key);
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
                log.warn("No default value for: " + key);
                return 0;
            }
            value = Float.valueOf(get(key, defValue));
            objectCache.put(key, value);
        }
        return value.floatValue();
    }


    /**
     * Get a property which is a delimited list of entries
     *
     * @param key
     * @return The string array of tokens, or an empty array if not present
     */
    private String[] getAsArray(String key) {
        String stringProp = get(key);
        if (stringProp == null) {
            return new String[0];
        } else {
            return stringProp.split(Globals.HISTORY_DELIMITER);
        }
    }

    public boolean hasExplicitValue(String key) {
        key = key.trim();
        return userPreferences.containsKey(key);
    }

    public void addOverrides(Map<String, String> newPrefs) {
        overrideKeys.addAll(newPrefs.keySet());
        userPreferences.putAll(newPrefs);
    }

    public boolean getAntiAliasing() {

        if (userPreferences.containsKey(Constants.ENABLE_ANTIALISING) || !Globals.IS_LINUX) {
            return getAsBoolean(Constants.ENABLE_ANTIALISING);
        } else {
            // Linux with no explicit setting
            return false;
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
            booleanCache.put(key, Boolean.valueOf(value));
        }
        objectCache.remove(key);
        mutationColorScheme = null;
    }

    private void clearCaches() {
        colorScaleCache.clear();
        booleanCache.clear();
        objectCache.clear();
        mutationColorScheme = null;
    }

    /**
     * Update preference.  This command is ignored if in batch mode *
     * @param key
     * @param value
     */
    public void put(String key, String value) {

        if (!Globals.isBatch()) {
            key = key.trim();

            // Explicitly setting removes override
            overrideKeys.remove(key);

            if (value == null || value.isBlank()) {
                userPreferences.remove(key);
            } else {
                userPreferences.put(key, value);
            }
            updateCaches(key, value);
            IGVEventBus.getInstance().post(new PreferencesChangeEvent());
        }
    }

    public void put(String key, boolean b) {
        put(key, String.valueOf(b));
    }


    public void putAll(Map<String, String> updatedPrefs) {
        for (Map.Entry<String, String> entry : updatedPrefs.entrySet()) {
            if (entry.getValue() == null || entry.getValue().isBlank()) {
                remove(entry.getKey());
            } else {
                put(entry.getKey(), entry.getValue());
            }
        }
        clearCaches();
        checkForProxyChanges(updatedPrefs);
        checkForAlignmentChanges(updatedPrefs);
        checkForCommandListenerChanges(updatedPrefs);
        checkForAttributePanelChanges(updatedPrefs);
        checkForCircViewChanges(updatedPrefs);
        checkForGoogleMenuChange(updatedPrefs);
        checkForRestartChanges(updatedPrefs);
        IGVEventBus.getInstance().post(new PreferencesChangeEvent());
    }

    private void checkForRestartChanges(Map<String, String> updatedPreferenceMap) {
        for(String key : RESTART_KEYS) {
            if(updatedPreferenceMap.containsKey(key)) {
                MessageUtils.showMessage("Preference changes will take effect after restart.");
                return;
            }
        }
    }

    private void checkForGoogleMenuChange(Map<String, String> updatedPreferenceMap) {

        if(updatedPreferenceMap.containsKey(ENABLE_GOOGLE_MENU) && IGV.hasInstance()) {
            try {
                IGVMenuBar.getInstance().enableGoogleMenu(getAsBoolean(ENABLE_GOOGLE_MENU));
            } catch (IOException e) {
                log.error("Error enabling/disabling Google menu", e);
            }
        }
    }

    private void checkForAlignmentChanges(Map<String, String> updatedPreferenceMap) {

        if (IGV.hasInstance()) {

            boolean reloadSAM = false;
            for (String key : SAM_RELOAD_KEYS) {
                if (updatedPreferenceMap.containsKey(key)) {
                    reloadSAM = true;
                    break;
                }
            }

            boolean refreshSAM = false;
            for (String key : SAM_REFRESH_KEYS) {
                if (updatedPreferenceMap.containsKey(key)) {
                    refreshSAM = true;
                    break;
                }
            }
            for (String key : BASEMOD_COLOR_KEYS) {
                if (updatedPreferenceMap.containsKey(key)) {
                    refreshSAM = true;
                    BaseModificationColors.updateColors();
                    break;
                }
            }

            for (String key : NUCLEOTIDE_COLOR_KEYS) {
                if (updatedPreferenceMap.containsKey(key)) {
                    SequenceRenderer.setNucleotideColors();
                    break;
                }
            }

            if (reloadSAM) {
                IGVEventBus.getInstance().post(new AlignmentTrackEvent(AlignmentTrackEvent.Type.RELOAD));
            }
            // A reload is harsher than a refresh; only send the weaker request if the stronger one is not sent.
            if (!reloadSAM && refreshSAM) {
                IGVEventBus.getInstance().post(new AlignmentTrackEvent(AlignmentTrackEvent.Type.REFRESH));
            }
            if (updatedPreferenceMap.containsKey(SAM_ALLELE_THRESHOLD)) {
                IGVEventBus.getInstance().post(new AlignmentTrackEvent(AlignmentTrackEvent.Type.ALLELE_THRESHOLD));
            }
        }

    }

    private void checkForProxyChanges(Map<String, String> updatedPreferenceMap) {
        if (HttpUtils.getInstance() != null) {
            for (String key : PROXY_KEYS) {
                if (updatedPreferenceMap.containsKey(key)) {
                    HttpUtils.getInstance().updateProxySettings();
                    break;
                }
            }
        }
    }

    private void checkForCommandListenerChanges(Map<String, String> updatedPreferenceMap) {
        if (updatedPreferenceMap.containsKey(PORT_ENABLED) || updatedPreferenceMap.containsKey(PORT_NUMBER)) {
            CommandListener.halt();
            if (getAsBoolean(PORT_ENABLED)) {
                CommandListener.start(getAsInt(PORT_NUMBER));
            }
        }
    }

    private void checkForAttributePanelChanges(Map<String, String> updatedPreferenceMap) {
        if (updatedPreferenceMap.containsKey(SHOW_ATTRIBUTE_VIEWS_KEY) || updatedPreferenceMap.containsKey(SHOW_DEFAULT_TRACK_ATTRIBUTES)) {
            if (IGV.hasInstance()) {
                IGV.getInstance().revalidateTrackPanels();
            }
        }
    }

    /**
     * Enabling circ view requires port listener
     */
    private void checkForCircViewChanges(Map<String, String> updatedPreferenceMap) {
        if (updatedPreferenceMap.containsKey(CIRC_VIEW_ENABLED) &&
                getAsBoolean(CIRC_VIEW_ENABLED) &&
                !getAsBoolean(PORT_ENABLED)) {
            put(PORT_ENABLED, true);
            CommandListener.start(getAsInt(PORT_NUMBER));
        }
    }


    public void remove(String key) {
        overrideKeys.remove(key);
        userPreferences.remove(key);
        booleanCache.remove(key);
        objectCache.remove(key);
        colorScaleCache.remove(key); //TODO same issue of cache not using String keys
        IGVEventBus.getInstance().post(new PreferencesChangeEvent());

    }


    public void clear() {
        userPreferences.clear();
        colorScaleCache.clear();
        booleanCache.clear();
        objectCache.clear();
        IGVEventBus.getInstance().post(new PreferencesChangeEvent());
    }


    public String getGenomeListURL() {
        return get(GENOMES_SERVER_URL);
    }

    public String getBackupGenomeListURL() {
        return get(BACKUP_GENOMES_SERVER_URL);
    }

    public String getProvisioningURL() {
        return get(PROVISIONING_URL);
    }

    public String getPortNumber() {
        return get(PORT_NUMBER);
    }

    public void overrideGenomeServerURL(String url) {
        userPreferences.put(GENOMES_SERVER_URL, url);
        overrideKeys.add(GENOMES_SERVER_URL);
        clearCaches();
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
    public void setLastSnapshotDirectory(String directory) {

        put(LAST_SNAPSHOT_DIRECTORY, directory);
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

        String lastFilePath = get(DEFINE_GENOME_INPUT_DIRECTORY_KEY, DirectoryManager.getUserDefaultDirectory().getAbsolutePath());

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

        String lastFilePath = get(LAST_GENOME_IMPORT_DIRECTORY, DirectoryManager.getUserDefaultDirectory().getAbsolutePath());

        if (lastFilePath != null) {
            genomeImportDirectory = new File(lastFilePath);
        }

        return genomeImportDirectory;
    }


    /**
     * @param bounds
     */
    public void setApplicationFrameBounds(Rectangle bounds) {

        if (bounds.width > 0 && bounds.height > 0) {
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
    }

    /**
     * @return
     */
    public Rectangle getApplicationFrameBounds() {

        Rectangle bounds = null;

        // Set the application's previous location and size
        String applicationBounds = get(FRAME_BOUNDS_KEY, null);

        if (applicationBounds != null) {
            String[] values = applicationBounds.split(",");
            int x = Integer.parseInt(values[0]);
            int y = Integer.parseInt(values[1]);
            int width = Integer.parseInt(values[2]);
            int height = Integer.parseInt(values[3]);

            if (width == 0 || height == 0) {
                return null;  // Don't know bounds
            }

            bounds = new Rectangle(x, y, width, height);
        }
        return bounds;
    }

    /**
     * @param recentSessions
     */
    public void setRecentSessions(String recentSessions) {
        remove(RECENT_SESSIONS);
        put(RECENT_SESSIONS, recentSessions);
    }


    public RecentFileSet getRecentSessions() {
        String sessionsString = get(RECENT_SESSIONS, null);
        return RecentFileSet.fromString(sessionsString, UIConstants.NUMBER_OF_RECENT_SESSIONS_TO_LIST);
    }

    public void setRecentUrls(String recentUrls) {
        remove(RECENT_URLS);
        put(RECENT_URLS, recentUrls);
    }


    public RecentUrlsSet getRecentUrls() {
        String sessionsString = get(RECENT_URLS, null);
        return RecentUrlsSet.fromString(sessionsString, UIConstants.NUMBER_OF_RECENT_SESSIONS_TO_LIST);
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
        userPreferences.put(DATA_SERVER_URL_KEY, url);
        overrideKeys.add(DATA_SERVER_URL_KEY);
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
        userPreferences.put(key, value);
        overrideKeys.add(key);    // <= order here is important, must be after userPreferences statement
        Map<String, String> updatedPrefs = new HashMap<String, String>();
        updatedPrefs.put(key, value);
        checkForAlignmentChanges(updatedPrefs);
        clearCaches();
    }

    public void setShowAttributeView(boolean isShowable) {
        put(SHOW_ATTRIBUTE_VIEWS_KEY, Boolean.toString(isShowable));
    }

    public void setLastGenome(String genomeId) {
        if (!genomeId.equals(get(DEFAULT_GENOME))) {
            put(DEFAULT_GENOME, genomeId);
        }
    }


    public String getDefaultGenome() {

        String genome = get(DEFAULT_GENOME, Globals.DEFAULT_GENOME);
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

    static String getMutationColorKey(String type) {
        switch (type) {
            case "Indel":
                return MUTATION_INDEL_COLOR_KEY;
            case "Missense":
                return MUTATION_MISSENSE_COLOR_KEY;
            case "Nonsense":
                return MUTATION_NONSENSE_COLOR_KEY;
            case "Splice_site":
                return MUTATION_SPLICE_SITE_COLOR_KEY;
            case "Synonymous":
                return MUTATION_SYNONYMOUS_COLOR_KEY;
            case "Targeted_region":
                return MUTATION_TARGETED_REGION_COLOR_KEY;
            case "Unknown":
                return MUTATION_UNKNOWN_COLOR_KEY;
            default:
                return "MUTATION_" + type + "_COLOR";
        }
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

    public static String generateGenomeIdString(Collection<GenomeListItem> genomeListItems) {
        String genomeString = "";

        for (GenomeListItem serverItem : genomeListItems) {
            genomeString += serverItem.getId() + Globals.HISTORY_DELIMITER;
        }

        genomeString = genomeString.substring(0, genomeString.length() - 1);
        return genomeString;
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
        return getAsArray(IGV_PLUGIN_LIST_KEY);
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

    public void print(PrintWriter pw) {
        for (Map.Entry<String, String> entry : userPreferences.entrySet()) {
            if (!overrideKeys.contains(entry.getKey())) {
                pw.print(entry.getKey());
                pw.print("=");
                pw.println(entry.getValue());
            }
        }
    }


}
