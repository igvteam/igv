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

package org.broad.igv.ui.color;

import org.broad.igv.Globals;
import org.broad.igv.session.Persistable;

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
public class PaletteColorTable implements ColorTable, Persistable {

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
