package org.broad.igv.ui.color;

import java.awt.*;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Experimental -- work in progress
 *
 * Table where Hue, Saturation, and Lightness are chosen randomly for each key
 *
 */
public class HSLColorTable extends ColorTable {

    final Map<String, Color> colorMap = new LinkedHashMap<>();

    final int hueCenter;
    boolean greyscale = false;


    public HSLColorTable(int hueCenter) {
        this.hueCenter = hueCenter;
    }

    protected Color computeColor(String key) {
        if (greyscale) {
            int r = (int) (50 + Math.random() * 155);
            return new Color(r, r, r);
        } else {
            double hue = ((hueCenter - 30) + Math.random() * 60);   // +/- 30 degrees
            double saturation = 0.4 + Math.random() * 0.6;     // [0.4, 1.0]
            double lightness = 0.2 + Math.random() * 0.4;      // [0.2, 0.6]
            int[] rgb = ColorUtilities.hslToRgb(hue, saturation, lightness);
            return new Color(rgb[0], rgb[1], rgb[2]);
        }
    }
}
