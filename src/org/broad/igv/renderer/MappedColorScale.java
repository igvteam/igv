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
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.renderer;

import org.broad.igv.ui.color.ColorUtilities;

import java.awt.*;
import java.util.HashMap;
import java.util.Map;

/**
 * A simple lookup color scale backed by a map.
 *
 * @author jrobinso
 */
public class MappedColorScale extends AbstractColorScale {

    public static String serializationClassId = "MappedColorScale";
    Map<String, Color> colorMap;
    boolean defaultCS = false;

    /**
     * Construct an instance from a serizled string representation.
     *
     * @param serializedInstance
     */
    public MappedColorScale(String serializedInstance) {

        colorMap = new HashMap<String, Color>();

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
        this.colorMap = new HashMap<String, Color>();
        this.colorMap.putAll(colorMap);
    }


    public boolean isDefault() {
        return defaultCS;
    }

    /**
     * Return a string representing the state fo this instance
     *
     * @return
     */
    public String asString() {
        StringBuffer buf = new StringBuffer();
        buf.append(serializationClassId);
        for (Map.Entry<String, Color> mapEntry : colorMap.entrySet()) {
            buf.append(";");
            buf.append(mapEntry.getKey() + " ");
            buf.append(ColorUtilities.colorToString(mapEntry.getValue()));

        }
        return buf.toString();
    }


    /**
     * @return the number of entries in this instance
     */
    public int getSize() {
        return colorMap.size();
    }

    /**
     * Comparison method.  Primary use is to support unit tests.
     *
     * @param anotherCS
     * @return
     */
    public boolean isSame(MappedColorScale anotherCS) {

        if (this.getSize() != anotherCS.getSize()) {
            return false;
        }
        for (String key : colorMap.keySet()) {
            if (!this.getColor(key).equals(anotherCS.getColor(key))) {
                return false;
            }
        }
        return true;

    }

    @Override
    public Color getColor(String key) {
        return colorMap.containsKey(key) ? colorMap.get(key) : defaultColor;
    }

    public Color getNoDataColor() {
        return noDataColor;
    }

}
