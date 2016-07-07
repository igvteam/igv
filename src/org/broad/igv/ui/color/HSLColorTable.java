package org.broad.igv.ui.color;

import java.awt.*;
import java.util.HashMap;
import java.util.Map;

/**
 * Experimental -- work in progress
 */
public class HSLColorTable implements ColorTable {

    int hueCenter;
    boolean greyscale = false;

    Map<String, Color> colorMap;

    public HSLColorTable(int hueCenter) {
        colorMap = new HashMap<>();
        this.hueCenter = hueCenter;
    }

    @Override
    public Color get(String key) {
        Color c = colorMap.get(key);
        if (c == null) {
            if (greyscale) {
                int r = (int) (50 + Math.random() * 155);
                c = new Color(r, r, r);
            } else {
                double hue = ((hueCenter - 30) + Math.random() * 60);   // +/- 30 degrees
                double saturation = 0.4 + Math.random() * 0.6;     // [0.4, 1.0]
                double lightness = 0.2 + Math.random() * 0.4;      // [0.2, 0.6]
                int[] rgb = ColorUtilities.hslToRgb(hue, saturation, lightness);
                c = new Color(rgb[0], rgb[1], rgb[2]);
            }

            colorMap.put(key, c);

        }
        return c;
    }
}
