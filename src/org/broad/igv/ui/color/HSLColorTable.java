package org.broad.igv.ui.color;

import sun.plugin2.util.ColorUtil;

import java.awt.*;
import java.util.HashMap;
import java.util.Map;

/**
 * Experimental -- work in progress
 */
public class HSLColorTable implements ColorTable {

//    int[] reds;
//    int[] greens;
//    int[] blues;
    int hueCenter;
    boolean greyscale = false;

    Map<String, Color> colorMap;


    public HSLColorTable(String type) {

        colorMap = new HashMap<>();

        if (type.equals("1")) {
            hueCenter = 30;

        } else if (type.equals("2")) {
            hueCenter = 270;
        } else {
            greyscale = true;
        }

    }

    @Override
    public Color get(String key) {
        Color c = colorMap.get(key);
        if (c == null) {
            if (greyscale) {
                int r =  (int) (50 + Math.random() * 155);
                c = new Color(r, r, r);
            } else {
                double hue = ((hueCenter - 30) + Math.random() * 60);   // +/- 30 degrees
                double saturation = 0.4 + Math.random() * 0.6;
                double lightness = 0.4 + Math.random() * 0.6;
                int[] rgb = ColorUtilities.hslToRgb(hue, saturation, lightness);
                c = new Color(rgb[0], rgb[1], rgb[2]);
            }

            colorMap.put(key, c);

        }
        return c;
    }
}
