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
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.*;
import org.broad.tribble.readers.AsciiLineReader;

import java.awt.*;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.io.File;
import java.io.IOException;
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
     * Map of data track identifiers (i.e. "array names") to its
     * attributeMap.   The attributeMap for a track maps attribute name (for
     * example "Cell Type"  to value (for example "ES");
     */
    LinkedHashMap<String, Map<String, String>> attributeMap = new LinkedHashMap();

    /**
     * List of attribute names (keys).  The list
     * is kept so the keys may be fetched in the order they were added.
     */
    List<String> attributeKeys = new ArrayList();

    /**
     * Column meta data (column == attributeKey).
     */
    Map<String, ColumnMetaData> columnMetaData = new HashMap();


    /**
     * The complete set of unique attribute values per attribute key.  This is useful in
     * assigning unique colors
     */
    Map<String, Set<String>> uniqueAttributeValues;

    Map<String, Color> colorMap = new Hashtable();

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
    public List<String> getAttributeKeys() {
        return new ArrayList(attributeKeys);
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
        final List<String> keys = getAttributeKeys();
        Set<String> hiddenKeys = IGV.getInstance().getSession().getHiddenAttributes();
        if (hiddenKeys != null) keys.removeAll(hiddenKeys);
        return keys;
    }

    public void clearAllAttributes() {
        attributeMap.clear();
        attributeKeys.clear();
        uniqueAttributeValues.clear();
        //hiddenAttributes.clear();
        loadedResources = new HashSet();
    }

    /**
     * Set the attribute value for the given track and key.
     */
    private void addAttribute(String trackIdentifier, String attributeName, String attributeValue) {

        if (attributeValue.equals("")) {
            return;
        }

        addAttributeKey(attributeName);

        Set<String> uniqueSet = uniqueAttributeValues.get(attributeName);
        if (uniqueSet == null) {
            uniqueSet = new HashSet<String>();
            uniqueAttributeValues.put(attributeName, uniqueSet);
        }
        uniqueSet.add(attributeValue);

        Map attributes = attributeMap.get(trackIdentifier);
        if (attributes == null) {
            attributes = new LinkedHashMap();
            attributeMap.put(trackIdentifier, attributes);
        }

        // attributeKey = column header, attributeValue = value for header
        // and track name (trackIdentifier) row intersection
        String key = attributeName.toUpperCase();
        attributes.put(key, attributeValue);
        updateMetaData(key, attributeValue);
    }

    public void addAttributeKey(String key) {
        if (!attributeKeys.contains(key) && !key.startsWith("#")) {
            attributeKeys.add(key);
        }
    }

    /**
     * Update the column meta data associated with the attribute key.
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

        if (metaData.isNumeric()) {
            try {
                double val = Double.parseDouble(attributeValue);
                metaData.updateRange(val);
            }
            catch (NumberFormatException e) {
                metaData.addNonNumericLabel(attributeValue);
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

        // If the file is "too large" better ask user
        // TODO -- ftp test
        final int oneMB = 1000000;
        long fileLength = 0;
        if (locator.isLocal()) {
            File f = new File(locator.getPath());
            fileLength = f.length();
        } else if (locator.getPath().startsWith("http")) {
            fileLength = IGVHttpClientUtils.getContentLength(new URL(locator.getPath()));
        }
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
            if (nextLine.toLowerCase().startsWith("#sampletable")) {
                loadSampleTable(reader, nextLine, locator.getPath());
            } else {
                loadOldSampleInfo(reader, nextLine, locator.getPath());
            }
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

    private void loadOldSampleInfo(AsciiLineReader reader, String nextLine, String path) throws IOException {
        // Parse column neadings for attribute names.
        // Columns 1 and 2 are array and sample name (not attributes)
        boolean foundAttributes = false;
        String[] colHeadings = nextLine.split("\t");
        int nLines = 0;
        int lineLimit = 100000;
        while ((nextLine = reader.readLine()) != null) {
            if (nLines++ > lineLimit) {
                break;
            }

            if (nextLine.startsWith("#colors")) {
                parseColors(reader);
                return;
            }

            String[] values = nextLine.split("\t");

            if (values.length >= 2) {
                String arrayName = values[0].trim();
                // Loop through attribute columns
                for (int i = 0; i < colHeadings.length; i++) {
                    String attributeName = colHeadings[i].trim();
                    String attributeValue = (i < values.length ? values[i].trim() : "");
                    addAttribute(arrayName, attributeName, attributeValue);
                    foundAttributes = true;
                }
            }
        }


        if (!foundAttributes) {
            throw new DataLoadException("Could not determine file type.  Does file have proper extension? ", path);
        }
    }

    private void parseColors(AsciiLineReader reader) throws IOException {

        String nextLine;
        while ((nextLine = reader.readLine()) != null) {
            String[] tokens = nextLine.split("\t");
            if (tokens.length >= 3) {
                String attKey = tokens[0];
                String attValue = tokens[1];
                Color color = ColorUtilities.stringToColor(tokens[2]);
                String key = (attKey + "_" + attValue).toLowerCase();
                colorMap.put(key, color);
            }
        }
    }

    /**
     * Load attributes from an ascii file in "Sample Table" format.  This format expects
     * a sample table section, prefaced by #sampleTable,  followed by a sample mapping
     * section  (track -> sample) prefaced by #sampleMappings
     */
    private void loadSampleTable(AsciiLineReader reader, String nextLine, String path) throws IOException {

        // Parse column neadings for attribute names.
        // Columns 1 and 2 are array and sample name (not attributes)
        nextLine = reader.readLine();
        String[] colHeadings = nextLine.split("\t");

        // Map of sample -> attribute list
        Map<String, List<Attribute>> sampleTable = new HashMap();

        boolean foundAttributes = false;
        int nLines = 0;
        int lineLimit = 100000;
        while ((nextLine = reader.readLine()) != null) {
            if (nLines++ > lineLimit || nextLine.toLowerCase().startsWith("#samplemapping")) {
                break;
            }
            String[] values = nextLine.split("\t");

            if (values.length >= 2) {
                String sampleName = values[0].trim();
                // Loop through attribute columns
                List<Attribute> attributes = new ArrayList(colHeadings.length);
                for (int i = 0; i < colHeadings.length; i++) {
                    String attributeName = colHeadings[i].trim();
                    String attributeValue = (i < values.length ? values[i].trim() : "");
                    attributes.add(new Attribute(attributeName, attributeValue));
                    foundAttributes = true;
                }
                sampleTable.put(sampleName, attributes);
            }
        }
        if (!foundAttributes) {
            throw new DataLoadException("Could not determine file type.  Does file have proper extension? ", path);
        }


        if (nextLine.toLowerCase().startsWith("#samplemapping")) {
            while ((nextLine = reader.readLine()) != null) {
                String[] tokens = nextLine.split("\t");
                if (tokens.length < 2) {
                    continue;
                }
                String array = tokens[0];
                String sample = tokens[1];
                List<Attribute> attributes = sampleTable.get(sample);
                if (attributes == null) {
                    log.info("Warning: sample in mapping section:  " + sample + " in sample table file " + path);
                } else {
                    for (Attribute att : attributes) {
                        addAttribute(array, att.getKey(), att.getValue());
                    }
                }
            }
        } else if (nextLine.startsWith("#colors")) {
            parseColors(reader);
            return;
        } else {
            // No mapping section.
            for (Map.Entry<String, List<Attribute>> entry : sampleTable.entrySet()) {
                String sample = entry.getKey();
                for (Attribute att : entry.getValue()) {
                    addAttribute(sample, att.getKey(), att.getValue());
                }
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

        String key = (attKey + "_" + attValue).toLowerCase();
        Color c = colorMap.get(key);
        if (c == null) {
            key = ("*_" + attValue).toLowerCase();
            c = colorMap.get(key);

            if (c == null) {
                Integer cnt = colorCounter.get(attKey);
                if (cnt == null) {
                    cnt = 0;
                }
                cnt++;
                colorCounter.put(attKey, cnt);
                float hue = (float) (.4 + 0.2 * Math.random());

                // int index = colorMap.size() + 1;
                c = randomColor(cnt);
                colorMap.put(key, c);
            }
        }
        return c;
    }


    static class ColumnMetaData {
        // Assume meta data is true until proven otherwise
        Set<String> nonNumericLabels = new HashSet();
        double min = Double.MAX_VALUE;
        double max = -min;

        void updateRange(double value) {
            min = Math.min(min, value);
            max = Math.max(max, value);
        }

        // Allow up to 1 non-numeric field

        public boolean isNumeric() {
            return nonNumericLabels.size() <= 1;
        }


        public void addNonNumericLabel(String label) {
            nonNumericLabels.add(label);

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
