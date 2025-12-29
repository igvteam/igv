/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.igv.renderer;

import org.igv.ui.color.ColorUtilities;

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
