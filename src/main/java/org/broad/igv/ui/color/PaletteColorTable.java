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
