/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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
 * AttributeManager.java
 *
 * Everything to do with attributes.
 */
package org.broad.igv.track;

import org.apache.log4j.Logger;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.renderer.AbstractColorScale;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.renderer.MonocolorScale;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.*;
import org.broad.tribble.readers.AsciiLineReader;

import java.awt.*;
import java.awt.font.NumericShaper;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.io.File;
import java.io.IOException;
import java.lang.reflect.Array;
import java.net.URL;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 */
public class AttributeManager {

    private static Logger log = Logger.getLogger(AttributeManager.class);

    private static AttributeManager singleton;
    final public static String ATTRIBUTES_LOADED_PROPERTY = "ATTRIBUTES_LOADED_PROPERTY";
    final public static String ATTRIBUTES_NARROWED_PROPERTY = "ATTRIBUTES_NARROWED_PROPERTY";

    private PropertyChangeSupport propertyChangeSupport;

    /**
     * The set of currently loaded attribute resource files
     */
    Set<ResourceLocator> loadedResources = new HashSet();


    /**
     * Sample table. Key is sample name, identifying a "row" in the table.  Value is a map of column name / value
     * pairs.   (e.g.  {TCGA-001  ->  { (gender->male),  (treated -> true), etc}}
     */
    LinkedHashMap<String, Map<String, String>> attributeMap = new LinkedHashMap();


    /**
     * Map of track id -> sample name.
     */
    Map<String, String> trackSampleMappings = new HashMap<String, String>();

    /**
     * List of attribute names.  The list
     * is kept so the keys may be fetched in the order they were added.
     */
    LinkedHashMap<String, String> attributeNames = new LinkedHashMap();

    /**
     * Column meta data (column == attributeKey).
     */
    Map<String, ColumnMetaData> columnMetaData = new HashMap();


    /**
     * The complete set of unique attribute values per attribute key.  This is useful in
     * assigning unique colors
     */
    Map<String, Set<String>> uniqueAttributeValues;

    /**
     * Maps symbolic (discrete) attribute values to colors. Key is a composite of attribute name and value
     */
    Map<String, Color> colorMap = new Hashtable();

    /**
     * Map of attribute column name -> color scale.   For numeric columns.
     */
    Map<String, AbstractColorScale> colorScales = new HashMap();

    Map<String, Integer> colorCounter = new HashMap();


    private AttributeManager() {
        propertyChangeSupport = new PropertyChangeSupport(this);
        uniqueAttributeValues = new HashMap();
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
        Map attributes = attributeMap.get(trackName);
        String key = attributeName.toUpperCase();
        return (attributes == null ? null : (String) attributes.get(key));
    }

    /**
     * Return the list of attribute names (keys) in the order they should
     * be displayed.
     */
    public List<String> getAttributeNames() {
        return new ArrayList(attributeNames.values());
    }

    /**
     * Return true if the associated column contains all numeric values
     */
    boolean isNumeric(String attributeName) {
        String key = attributeName.toUpperCase();
        ColumnMetaData metaData = columnMetaData.get(key);
        return metaData != null && metaData.isNumeric();
    }


    // TODO -- don't compute this on the fly every time its called

    public List<String> getVisibleAttributes() {
        final Set<String> allKeys = attributeNames.keySet();
        Set<String> hiddenAttributes = IGV.getInstance().getSession().getHiddenAttributes();
        if (hiddenAttributes != null) {
            allKeys.removeAll(hiddenAttributes);
        }

        ArrayList<String> visibleAttributes = new ArrayList<String>(allKeys.size());
        for (String key : allKeys) {
            visibleAttributes.add(attributeNames.get(key));
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
     * Set the attribute value for the given track or sample id and key.
     */
    public void addAttribute(String trackIdentifier, String name, String attributeValue) {

        if (attributeValue.equals("")) {
            return;
        }

        String key = name.toUpperCase();
        addAttributeName(name);

        Set<String> uniqueSet = uniqueAttributeValues.get(key);
        if (uniqueSet == null) {
            uniqueSet = new HashSet<String>();
            uniqueAttributeValues.put(key, uniqueSet);
        }
        uniqueSet.add(attributeValue);

        Map attributes = attributeMap.get(trackIdentifier);
        if (attributes == null) {
            attributes = new LinkedHashMap();
            attributeMap.put(trackIdentifier, attributes);
        }

        // attributeKey = column header, attributeValue = value for header
        // and track name (trackIdentifier) row intersection
        attributes.put(key, attributeValue);
        updateMetaData(key, attributeValue);
    }

    public void addAttributeName(String name) {
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
            metaData = new ColumnMetaData();
            columnMetaData.put(key, metaData);
        }

        // Test if data is numeric.  Skip null and blank values
        if (attributeValue != null && attributeValue.length() > 0 && metaData.isNumeric()) {
            try {
                double val = Double.parseDouble(attributeValue);
                metaData.updateRange(val);
            } catch (NumberFormatException e) {
                metaData.markNonNumeric();
            }
        }


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

        // If the file is "too large"  ask user
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
        String nextLine = null;
        try {
            reader = ParsingUtils.openAsciiReader(locator);
            nextLine = reader.readLine();
            loadSampleTable(reader, nextLine, locator.getPath());

            loadedResources.add(locator);

            //createCurrentAttributeFileString(files);
            IGV.getInstance().getTrackManager().resetOverlayTracks();

            IGV.getInstance().doRefresh();

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

    /**
     * Load sample table, which might optionally have 3 sections
     * #sampletable (default)
     * #samplemappint (track id -> sample mapping table)
     * #colors (color table)
     */
    private void loadSampleTable(AsciiLineReader reader, String nextLine, String path) throws IOException {

        String[] colHeadings = null;
        
        List<String> sections = Arrays.asList("#sampletable", "#samplemapping", "#colors");

        boolean foundAttributes = false;
        int nLines = 0;
        int lineLimit = 100000;
        String section = "#sampletable";
        while ((nextLine = reader.readLine()) != null) {
            if (nLines++ > lineLimit) {
                break;
            }
            if (nextLine.toLowerCase().startsWith("#")) {
                String tmp =  nextLine.toLowerCase().trim();
                if(sections.contains(tmp)) {
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
        return trackSampleMappings.get(track);
    }

    /**
     * Represents a specific attribute instance.
     *
     * @author jrobinso
     */
    private static class Attribute {
        private String key;
        private String value;

        public Attribute(String key, String value) {
            this.key = key;
            this.value = value;
        }

        public String getKey() {
            return key;
        }

        public String getValue() {
            return value;
        }

    }


    public Color getColor(String attKey, String attValue) {

        if (attValue == null || attValue.length() == 0) {
            return Color.white;
        }

        if (isNumeric(attKey)) {
            AbstractColorScale cs = colorScales.get(attKey);
            {
                if (cs != null) {
                    try {
                        float x = Float.parseFloat(attValue);
                        return cs.getColor(x);
                    } catch (NumberFormatException e) {
                        return Color.white;
                    }
                }

            }
        }

        String key = (attKey + "_" + attValue).toUpperCase();
        Color c = colorMap.get(key);
        if (c == null) {
            key = ("*_" + attValue).toUpperCase();
            c = colorMap.get(key);

            if (c == null) {

                key = (attValue + "_*").toUpperCase();
                c = colorMap.get(key);

                if (c == null) {

                    Integer cnt = colorCounter.get(attKey);
                    if (cnt == null) {
                        cnt = 0;
                    }
                    cnt++;
                    colorCounter.put(attKey, cnt);
                    c = randomColor(cnt);
                }

                colorMap.put(key, c);
            }
        }
        return c;
    }


    static class ColumnMetaData {
        // Assume meta data is true until proven otherwise
        boolean numeric = true;
        double min = Double.MAX_VALUE;
        double max = -min;

        void updateRange(double value) {
            min = Math.min(min, value);
            max = Math.max(max, value);
        }

        // Allow up to 1 non-numeric field

        public boolean isNumeric() {
            return numeric;
        }


        public void markNonNumeric() {
            numeric = false;
        }
    }

    static class Range {
        double min;
        double max;
        Color color;
    }


    public static Color randomColor(int idx) {
        float hue = (float) Math.random();
        float sat = (float) (0.8 * Math.random());
        float bri = (float) (0.6 + 0.4 * Math.random());
        return Color.getHSBColor(hue, sat, bri);
    }

}
