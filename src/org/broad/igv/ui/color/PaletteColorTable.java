/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.ui.color;

import org.broad.igv.Globals;

import java.awt.*;
import java.util.*;

/**
 * A lookup table mapping symbols (strings) -> color.  Can be initiated with our without a palette.  If
 * a palette is used colors are selected sequentially from the palette as needed until it is exhausted.
 * <p/>
 *
 * @author jrobinso
 * @date Feb 26, 2010
 */
public class PaletteColorTable implements ColorTable {

    LinkedHashMap<String, Color> colorMap;
    Color[] colors;

    public PaletteColorTable() {
        colorMap = new LinkedHashMap();
    }

    public PaletteColorTable(ColorPalette palette) {
        if (palette != null) {
            this.colors = palette.getColors();
        }
        colorMap = new LinkedHashMap();
    }

    public void put(String key, Color c) {
        colorMap.put(key, c);
    }

    public Color get(String key) {
        Color c = colorMap.get(key);
        if (c == null) {
            final int colorIdx = colorMap.size();
            if (colors != null && colorIdx < colors.length) {
                c = colors[colorIdx];
            } else {
                c = ColorUtilities.randomColor(colorIdx);
            }
            colorMap.put(key, c);
        }
        return c;
    }

    public Collection<String> getKeys() {
        return colorMap.keySet();
    }

    public Set<Map.Entry<String, Color>> entrySet() {
        return colorMap.entrySet();
    }


    public String getMapAsString() {
        StringBuffer buf = new StringBuffer();
        boolean firstEntry = true;
        for (Map.Entry<String, Color> entry : colorMap.entrySet()) {
            if (!firstEntry) {
                buf.append(";");
            }
            String cs = ColorUtilities.colorToString(entry.getValue());
            buf.append(entry.getKey());
            buf.append("=");
            buf.append(cs);
            firstEntry = false;
        }
        return buf.toString();
    }

    public void restoreMapFromString(String string) {

        if (string == null || string.isEmpty()) return;
        colorMap.clear();
        String[] tokens = Globals.semicolonPattern.split(string);
        for (String t : tokens) {
            String[] kv = Globals.equalPattern.split(t);
            colorMap.put(kv[0], ColorUtilities.stringToColor(kv[1]));
        }
    }


    /**
     * Return object state as a map of key-value string pairs
     *
     * @return
     */
    public Map<String, String> getPersistentState() {
        Map<String, String> state = new HashMap<String, String>();
        if (colorMap != null && colorMap.size() > 0) {
            state.put("colorMap", getMapAsString());
        }
        return state;
    }

    /**
     * Restore object state from a map of key-value string pairs
     *
     * @param values
     */
    public void restorePersistentState(Map<String, String> values) {
        String colorMapString = values.get("colorMap");
        if (colorMapString != null) {
            restoreMapFromString(colorMapString);
        }
    }

    public LinkedHashMap<String,Color> getColorMap() {
        return colorMap;
    }
}
