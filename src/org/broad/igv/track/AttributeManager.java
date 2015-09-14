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
 * AttributeManager.java
 *
 * Everything to do with attributes.
 */
package org.broad.igv.track;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.renderer.AbstractColorScale;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.renderer.MonocolorScale;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.color.ColorPalette;
import org.broad.igv.ui.color.ColorTable;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.color.PaletteColorTable;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.Utilities;
import htsjdk.tribble.readers.AsciiLineReader;

import java.awt.*;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.io.IOException;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 */
public class AttributeManager {

    public static List<String> defaultTrackAttributes = Arrays.asList(Globals.TRACK_NAME_ATTRIBUTE,
            Globals.TRACK_DATA_FILE_ATTRIBUTE, Globals.TRACK_DATA_TYPE_ATTRIBUTE);
    private static Logger log = Logger.getLogger(AttributeManager.class);

    private static AttributeManager singleton;
    final public static String ATTRIBUTES_LOADED_PROPERTY = "ATTRIBUTES_LOADED_PROPERTY";
    final public static String ATTRIBUTES_NARROWED_PROPERTY = "ATTRIBUTES_NARROWED_PROPERTY";

    private PropertyChangeSupport propertyChangeSupport;

    /**
     * The set of currently loaded attribute resource files
     */
    Set<ResourceLocator> loadedResources = Collections.synchronizedSet(new HashSet());


    /**
     * Sample table. Key is sample name, identifying a "row" in the table.  Value is a map of column name / value
     * pairs.   (e.g.  {TCGA-001  ->  { (gender->male),  (treated -> true), etc}}
     */
    Map<String, Map<String, String>> attributeMap = Collections.synchronizedMap(new LinkedHashMap());


    /**
     * Map of track id -> sample name.
     */
    Map<String, String> trackSampleMappings = Collections.synchronizedMap(new HashMap<String, String>());

    /**
     * List of attribute names.  The list
     * is kept so the keys may be fetched in the order they were added.
     */
    Map<String, String> attributeNames = Collections.synchronizedMap(new LinkedHashMap());

    /**
     * Column meta data (column == attributeKey).
     */
    Map<String, ColumnMetaData> columnMetaData = Collections.synchronizedMap(new HashMap());


    /**
     * The complete set of unique attribute values per attribute key.  This is useful in
     * assigning unique colors
     */
    Map<String, Set<String>> uniqueAttributeValues = Collections.synchronizedMap(new HashMap());

    /**
     * Maps symbolic (discrete) attribute values to colors. Key is a composite of attribute name and value
     */
    Map<String, Color> colorMap = Collections.synchronizedMap(new HashMap());

    /**
     * Map of attribute column name -> color scale.   For numeric columns.
     */
    Map<String, AbstractColorScale> colorScales = new HashMap();

    Map<String, ColorTable> colorTables = new HashMap<String, ColorTable>();

    Map<String, Integer> colorCounter = new HashMap();

    private AttributeManager() {
        propertyChangeSupport = new PropertyChangeSupport(this);
        //hiddenAttributes.add("NAME");
        //hiddenAttributes.add("DATA FILE");
        //hiddenAttributes.add("DATA TYPE");

    }

    static synchronized public AttributeManager getInstance() {

        if (singleton == null) {
            singleton = new AttributeManager();
        }
        return singleton;
    }

    public void addPropertyChangeListener(PropertyChangeListener listener) {
        propertyChangeSupport.addPropertyChangeListener(listener);
    }

    public void removePropertyChangeListener(PropertyChangeListener listener) {
        propertyChangeSupport.removePropertyChangeListener(listener);
    }

    /**
     * Return the attribute value for the given track (trackName) and key.
     */
    public String getAttribute(String trackName, String attributeName) {
        Map<String, String> attributes = attributeMap.get(trackName);
        String key = attributeName.toUpperCase();
        String value = attributes == null ? null : attributes.get(key);
        if (value == null && trackSampleMappings.containsKey(trackName)) {
            final String sample = trackSampleMappings.get(trackName);
            attributes = attributeMap.get(sample);
            if (attributes != null) {
                value = attributes.get(key);
            }
        }
        return value;
    }

    /**
     * Return the list of attribute names (keys) in the order they should
     * be displayed.
     */
    public List<String> getAttributeNames() {
        ArrayList<String> attNames = new ArrayList<String>(attributeNames.values());
        return attNames;
    }

    /**
     * Return true if the associated column contains all numeric values
     */
    public boolean isNumeric(String attributeName) {
        String key = attributeName.toUpperCase();
        ColumnMetaData metaData = columnMetaData.get(key);
        return metaData != null && metaData.isNumeric();
    }


    /**
     * Return all attributes, except those that have been "hidden" in the attribute panel
     * TODO -- don't compute this every time (or at least profile to see if this is a problem).
     *
     * @return
     */
    public List<String> getVisibleAttributes() {

        List<String> visibleAttributes = getAttributeNames();
        if (visibleAttributes == null) {
            Collections.emptyList();
        }
        final Set<String> hiddenAttributes = IGV.getInstance().getSession().getHiddenAttributes();
        if (hiddenAttributes != null) {
            visibleAttributes.removeAll(hiddenAttributes);
        }

        return visibleAttributes;
    }

    public void clearAllAttributes() {
        attributeMap.clear();
        attributeNames.clear();
        uniqueAttributeValues.clear();
        //hiddenAttributes.clear();
        loadedResources = new HashSet();
    }

    /**
     * Set an attribute value
     *
     * @param rowId          -- track or sample identifier
     * @param attributeName
     * @param attributeValue
     */
    public void addAttribute(String rowId, String attributeName, String attributeValue) {

        if (attributeValue == null || attributeValue.equals("")) {
            return;
        }

        // Add the 3 "special" attributes to ensure they are the first columns
        if (attributeNames.isEmpty()) {
            addAttributeName("NAME");
            addAttributeName("DATA TYPE");
            addAttributeName("DATA FILE");
        }

        addAttributeName(attributeName);

        String key = attributeName.toUpperCase();
        Set<String> uniqueSet = uniqueAttributeValues.get(key);
        if (uniqueSet == null) {
            uniqueSet = new HashSet<String>();
            uniqueAttributeValues.put(key, uniqueSet);
        }
        uniqueSet.add(attributeValue);

        Map attributes = attributeMap.get(rowId);
        if (attributes == null) {
            attributes = new LinkedHashMap();
            attributeMap.put(rowId, attributes);
        }

        // attributeKey = column header, attributeValue = value for header
        // and track name (trackIdentifier) row intersection
        attributes.put(key, attributeValue);
        updateMetaData(key, attributeValue);
    }

    private void addAttributeName(String name) {
        String key = name.toUpperCase();
        if (!attributeNames.containsKey(key) && !name.startsWith("#")) {
            attributeNames.put(key, name);
        }
    }


    /**
     * Update the column meta data associated with the attribute key.
     * <p/>
     * Note: Currently the meta data only records if the column is numeric.
     *
     * @param attributeName
     * @param attributeValue
     */
    private void updateMetaData(String attributeName, String attributeValue) {

        String key = attributeName.toUpperCase();
        ColumnMetaData metaData = columnMetaData.get(key);
        if (metaData == null) {
            metaData = new ColumnMetaData(key);
            columnMetaData.put(key, metaData);
        }

        metaData.updateMetrics(attributeValue);


    }

    /**
     * Test to see if this file could be a sample information file.  Some characteristics are (1) is tab delimited
     * with at least 2 columns,  (2) is ascii,  (3) is not too large
     *
     * @param locator
     * @return
     */
    public static boolean isSampleInfoFile(ResourceLocator locator) throws IOException {

        if (!FileUtils.isTabDelimited(locator, 2)) {
            return false;
        }

        // If the file is too large, ask user
        // TODO -- ftp test
        final int oneMB = 1000000;
        long fileLength = ParsingUtils.getContentLength(locator.getPath());
        if (fileLength > oneMB) {
            return MessageUtils.confirm("<html>Cannot determine file type of: " + locator.getPath() +
                    "<br>Is this a sample information file?");
        }
        return true;
    }

    /**
     * Load attributes from an ascii file in "Sample Info" format.
     */
    public void loadSampleInfo(ResourceLocator locator) {
        AsciiLineReader reader = null;
        try {
            reader = ParsingUtils.openAsciiReader(locator);

            loadSampleTable(reader, locator.getPath());

            loadedResources.add(locator);

            if (!Globals.isHeadless()) {
                IGV.getInstance().resetOverlayTracks();
                IGV.getInstance().doRefresh();
            }

        } catch (IOException ex) {
            log.error("Error loading attribute file", ex);
            throw new DataLoadException("Error reading attribute file", locator.getPath());
        } finally {
            if (reader != null) {
                reader.close();

            }
            firePropertyChange(this, ATTRIBUTES_LOADED_PROPERTY, null, null);
        }
    }


    static Set<String> nonGroupable = new HashSet<String>(Arrays.asList("DATA FILE", "DATA TYPE",
            "VITALSTATUS", "VITAL STATUS", "KARNSCORE", "CENSURED"));

    public List<String> getGroupableAttributes() {
        List<String> seriesNames = new ArrayList<String>();
        for (Map.Entry<String, Set<String>> entry : uniqueAttributeValues.entrySet()) {
            int cnt = entry.getValue().size();
            String att = entry.getKey();
            if (cnt > 1 && cnt < 10 && !nonGroupable.contains(att)) {
                seriesNames.add(att);
            }
        }

        return seriesNames;

    }

    /**
     * Load sample table, which might optionally have 3 sections
     * #sampletable (default)
     * #samplemappint (track id -> sample mapping table)
     * #colors (color table)
     */
    private void loadSampleTable(AsciiLineReader reader, String path) throws IOException {

        String[] colHeadings = null;

        List<String> sections = Arrays.asList("#sampletable", "#samplemapping", "#colors");

        boolean foundAttributes = false;
        int nLines = 0;
        int lineLimit = 100000;
        String nextLine;
        String section = "#sampletable";
        while ((nextLine = reader.readLine()) != null) {
            if (nLines++ > lineLimit) {
                break;
            }
            if (nextLine.toLowerCase().startsWith("#")) {
                String tmp = nextLine.toLowerCase().trim();
                if (sections.contains(tmp)) {
                    section = tmp;
                }
                continue;
            }

            String[] tokens = nextLine.split("\t");

            if (section.equals("#sampletable")) {
                if (tokens.length >= 2) {
                    if (colHeadings == null) {
                        colHeadings = tokens;
                    } else {
                        String sampleName = tokens[0].trim();

                        // Loop through attribute columns
                        //List<Attribute> attributes = new ArrayList(colHeadings.length);
                        for (int i = 0; i < colHeadings.length; i++) {
                            String attributeName = colHeadings[i].trim();
                            String attributeValue = (i < tokens.length ? tokens[i].trim() : "");
                            addAttribute(sampleName, attributeName, attributeValue);
                            foundAttributes = true;
                        }
                    }
                }
            } else if (section.equals("#samplemapping")) {
                if (tokens.length < 2) {
                    continue;
                }
                String track = tokens[0];
                String sample = tokens[1];
                trackSampleMappings.put(track, sample);

            } else if (section.equals("#colors")) {
                parseColors(tokens);

            }
        }
        if (!foundAttributes) {
            throw new DataLoadException("Could not determine file type.  Does file have proper extension? ", path);
        }


    }


    private void parseColors(String[] tokens) throws IOException {

        if (tokens.length >= 3) {
            String attKey = tokens[0].toUpperCase();
            if (isNumeric(attKey)) {

                ColumnMetaData metaData = columnMetaData.get(attKey);
                String rangeString = tokens[1].trim();
                float min = (float) metaData.min;
                float max = (float) metaData.max;
                if (!rangeString.equals("*") && rangeString.length() > 0) {
                    String[] tmp = rangeString.split(":");
                    if (tmp.length > 1) {
                        try {
                            min = Float.parseFloat(tmp[0]);
                            max = Float.parseFloat(tmp[1]);
                        } catch (NumberFormatException e) {
                            log.error("Error parsing range string: " + rangeString, e);
                        }
                    }
                }

                AbstractColorScale scale = null;
                if (tokens.length == 3) {
                    Color baseColor = ColorUtilities.stringToColor(tokens[2]);
                    scale = new MonocolorScale(min, max, baseColor);
                    colorScales.put(attKey, scale);
                } else {
                    Color color1 = ColorUtilities.stringToColor(tokens[2]);
                    Color color2 = ColorUtilities.stringToColor(tokens[3]);
                    if (min < 0) {
                        scale = new ContinuousColorScale(min, 0, max, color1, Color.white, color2);
                    } else {
                        scale = new ContinuousColorScale(min, max, color1, color2);
                    }
                }
                colorScales.put(attKey, scale);

            } else {
                String attValue = tokens[1];
                Color color = ColorUtilities.stringToColor(tokens[2]);
                String key = (attKey + "_" + attValue).toUpperCase();
                colorMap.put(key, color);
            }

        }
    }


    public void firePropertyChange(Object source, String propertyName,
                                   Object oldValue, Object newValue) {

        PropertyChangeEvent event =
                new PropertyChangeEvent(
                        source,
                        propertyName,
                        oldValue,
                        newValue);
        propertyChangeSupport.firePropertyChange(event);
    }

    public Comparator getAttributeComparator() {
        return Utilities.getNumericStringComparator();
    }

    /**
     * @return set of curently loaded resources
     */
    public Set<ResourceLocator> getLoadedResources() {
        return loadedResources;
    }

    public String getSampleFor(String track) {

        //return trackSampleMappings.get(track);

        if (trackSampleMappings.containsKey(track)) {
            return trackSampleMappings.get(track);
        } else if (isTCGAName(track)) {
            String sample = track.substring(0, 12);
            addAttribute(track, "Sample", sample);
            trackSampleMappings.put(track, sample);
            return sample;
        } else {
            String key = PreferenceManager.getInstance().get(PreferenceManager.OVERLAY_ATTRIBUTE_KEY);
            return key == null ? null : getAttribute(track, key);
        }
    }

    // TCGA identifers have the form TCGA-00-0000
    public static boolean isTCGAName(String name) {
        return name.length() >= 12 && name.toUpperCase().startsWith("TCGA-") &&
                name.charAt(7) == '-';

    }

    public Color getColor(String attKey, String attValue) {

        if (attValue == null || attValue.length() == 0) {
            return Color.gray;
        }

        final ColumnMetaData metaData = columnMetaData.get(attKey.toUpperCase());
        if (metaData == null) {
            return Color.gray;
        }
        if (metaData.isNumeric()) {
            AbstractColorScale cs = colorScales.get(attKey);
            {
                if (cs == null) {
                    // Create color scale based loosely on Brewer diverging / sequential palletes
                    // TODO -- use actual brewer palletes if # of values < 8
                    if (metaData.isDiverging()) {
                        // reg-blue diverging
                        Color minColor = new Color(198, 219, 239);
                        Color midColor = Color.white;
                        Color maxColor = new Color(33, 102, 172);
                        cs = new ContinuousColorScale(metaData.getMin(), 0, metaData.getMax(), minColor, midColor, maxColor);
                        colorScales.put(attKey, cs);

                    } else {
                        // Blues scale
                        Color minColor = new Color(198, 219, 239);
                        Color maxColor = new Color(8, 69, 148);
                        cs = new ContinuousColorScale(metaData.getMin(), metaData.getMax(), minColor, maxColor);
                        colorScales.put(attKey, cs);
                    }
                }
                try {
                    float x = Float.parseFloat(attValue);
                    return cs.getColor(x);
                } catch (NumberFormatException e) {
                    return Color.lightGray;
                }


            }
        }

        // Look for color in pre-loaded color map
        String key = (attKey + "_" + attValue).toUpperCase();
        Color c = colorMap.get(key);
        if (c == null) {
            key = ("*_" + attValue).toUpperCase();
            c = colorMap.get(key);
            if (c == null) {
                key = (attValue + "_*").toUpperCase();
                c = colorMap.get(key);
            }
        }

        // Get color from palette
        if (c == null) {

            // Measure of "information content" added by using color, very crude
            //boolean useColor = (metaData.getUniqueCount() < 10 || metaData.getUniqueRatio() <= 0.5) &&
            //        !(attKey.equals("NAME") || attKey.equals("DATA FILE") || attKey.equals("DATA TYPE"));
            boolean useColor = true;
            if (useColor) {
                ColorTable ct = colorTables.get(attKey);
                if (ct == null) {
                    ColorPalette palette = ColorUtilities.getNextPalette();
                    ct = new PaletteColorTable(palette);
                    colorTables.put(attKey, ct);
                }
                c = ct.get(attValue);
            } else {
                c = ColorUtilities.randomDesaturatedColor(0.5f);
                colorMap.put(key, c);
            }
        }

        return c;
    }


    public ColumnMetaData getColumnMetaData(String key) {
        return columnMetaData.get(key.toUpperCase());
    }

    public static class ColumnMetaData {

        String name;
        private double min = Double.MAX_VALUE;
        private double max = -min;

        int totalCount = 0;
        public HashSet<String> uniqueAlphaValues = new HashSet<String>();
        public HashSet<String> uniqueNumericValues = new HashSet<String>();

        ColumnMetaData(String name) {
            this.name = name;
        }

        public void updateMetrics(String attributeValue) {

            totalCount++;

            // Test if data is numeric.  Skip null and blank values
            if (attributeValue != null && attributeValue.length() > 0) {
                try {
                    double value = Double.parseDouble(attributeValue);
                    uniqueNumericValues.add(attributeValue);
                    min = Math.min(min, value);
                    max = Math.max(max, value);

                } catch (NumberFormatException e) {
                    uniqueAlphaValues.add(attributeValue);
                }
            }
        }

        /**
         * A column is considered numeric if it has at least 2 numeric values, and
         * no more than 1 non-numeric value.
         *
         * @return
         */
        public boolean isNumeric() {
            return uniqueNumericValues.size() > 1 && uniqueAlphaValues.size() < 2;
        }

        public boolean isDiverging() {
            return min < 0;
        }

        public double getMin() {
            return min;
        }

        public double getMax() {
            return max;
        }

        public double getUniqueRatio() {

            double totalUnique = uniqueAlphaValues.size() + uniqueNumericValues.size();
            return totalUnique / totalCount;

        }

        public int getUniqueCount() {
            return uniqueAlphaValues.size() + uniqueNumericValues.size();
        }

        public int getTotalCount() {
            return totalCount;
        }
    }


}
