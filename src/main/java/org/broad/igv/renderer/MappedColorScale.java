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

package org.broad.igv.renderer;

import org.broad.igv.ui.color.ColorUtilities;

import java.awt.*;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * A simple lookup color scale backed by a map.
 *
 * @author jrobinso
 */
public class MappedColorScale extends AbstractColorScale {

    public static final String serializationClassId = "MappedColorScale";
    final Map<String, Color> colorMap = new LinkedHashMap<>();

    /**
     * Construct an instance from a serialized string representation.
     *
     * @param serializedInstance
     */
    public MappedColorScale(String serializedInstance) {

        String[] tokens = serializedInstance.split(";");

        // First token is class string -- todo, validate this?
        for (int i = 1; i < tokens.length; i++) {
            String[] keyValue = tokens[i].split(" ");
            if (keyValue.length == 2) {
                String key = keyValue[0].trim();
                Color c = ColorUtilities.stringToColor(keyValue[1]);
                colorMap.put(key, c);
            }
        }
    }

    /**
     * Construct an instance from a color map
     *
     * @param colorMap
     */
    public MappedColorScale(Map<String, Color> colorMap) {
        this.colorMap.putAll(colorMap);
    }


    /**
     * Return a string representing the state of this instance
     *
     * @return
     */
    public String asString() {
        StringBuilder buf = new StringBuilder();
        buf.append(serializationClassId);
        colorMap.forEach((key, value) -> {
            buf.append(";");
            buf.append(key + " ");
            buf.append(ColorUtilities.colorToString(value));
        });
        return buf.toString();
    }


    /**
     * @return the number of entries in this instance
     */
    public int getSize() {
        return colorMap.size();
    }

    @Override
    public Color getColor(String key) {
        return colorMap.getOrDefault(key, DEFAULT_COLOR);
    }

}
