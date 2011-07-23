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


package org.broad.igv.util;


import org.apache.log4j.Logger;

import java.awt.*;
import java.util.HashMap;
import java.util.Map;
import java.util.WeakHashMap;


/**
 * Miscellaneous utilities for parsing and manipulating colors.
 *
 * @author Jim Robinson
 */
public class ColorUtilities {

    private static Logger log = Logger.getLogger(ColorUtilities.class);

    public static Map<Object, Color> colorCache = new WeakHashMap<Object, Color>(100);

    private static float[] whiteComponents = Color.white.getRGBColorComponents(null);

    // HTML 4.1 color table,  + orange and magenta
    static Map<String, String> colorSymbols = new HashMap();
    static {
        colorSymbols.put("white", "FFFFFF");
        colorSymbols.put("silver", "C0C0C0");
        colorSymbols.put("gray", "808080");
        colorSymbols.put("black", "000000");
        colorSymbols.put("red", "FF0000");
        colorSymbols.put("maroon", "800000");
        colorSymbols.put("yellow", "FFFF00");
        colorSymbols.put("olive", "808000");
        colorSymbols.put("lime", "00FF00");
        colorSymbols.put("green", "008000");
        colorSymbols.put("aqua", "00FFFF");
        colorSymbols.put("teal", "008080");
        colorSymbols.put("blue", "0000FF");
        colorSymbols.put("navy", "000080");
        colorSymbols.put("fuchsia", "FF00FF");
        colorSymbols.put("purple", "800080");
        colorSymbols.put("orange", "FFA500");
        colorSymbols.put("magenta", "FF00FF");
    }


    /**
     * Method description
     *
     * @param idx
     * @param alpha
     * @return
     */
    public static Color randomColor(int idx, float alpha) {

        int col1 = 0;
        int col2 = 0;
        int col3 = 0;

        int BASE_COL = 40;
        int RAND_COL = 255 - BASE_COL;

        idx += 1;    // avoid 0
        col1 = Math.abs(BASE_COL + (idx * 33) % RAND_COL);
        col2 = Math.abs(BASE_COL + (idx * 55) % RAND_COL);
        col3 = Math.abs(BASE_COL + (idx * 77) % RAND_COL);

        return new Color(col1, col2, col3, (int) (255 * alpha));
    }

    /**
     * Port of DChip function of the same name.
     *
     * @param idx
     * @return
     */
    public static Color randomColor(int idx) {

        int col1 = 0;
        int col2 = 0;
        int col3 = 0;

        int BASE_COL = 40;
        int RAND_COL = 255 - BASE_COL;

        idx += 1;    // avoid 0
        col1 = Math.abs(BASE_COL + (idx * 33) % RAND_COL);
        col2 = Math.abs(BASE_COL + (idx * 55) % RAND_COL);
        col3 = Math.abs(BASE_COL + (idx * 77) % RAND_COL);

        return new Color(col1, col2, col3);
    }


    private static float[] hsbvals = new float[3];

    /**
     * Method description
     *
     * @param inputColor
     * @param hue
     * @param saturation
     * @param brightness
     * @return
     */
    public static Color adjustHSB(Color inputColor, float hue, float saturation, float brightness) {
        Color.RGBtoHSB(inputColor.getRed(), inputColor.getGreen(), inputColor.getBlue(), hsbvals);
        return Color.getHSBColor(hue * hsbvals[0], saturation * hsbvals[1],
                brightness * hsbvals[2]);
    }


    public static String colorToString(Color color) {

        StringBuffer buffer = new StringBuffer();
        buffer.append(color.getRed());
        buffer.append(",");
        buffer.append(color.getGreen());
        buffer.append(",");
        buffer.append(color.getBlue());
        return buffer.toString();
    }

    public static Color stringToColor(String string) {

        try {
            Color c = colorCache.get(string);
            if (c == null) {
                if (string.contains(",")) {
                    String[] rgb = string.split(",");
                    int red = Integer.parseInt(rgb[0]);
                    int green = Integer.parseInt(rgb[1]);
                    int blue = Integer.parseInt(rgb[2]);
                    c = new Color(red, green, blue);
                } else if (string.startsWith("#")) {
                    c = hexToColor(string.substring(1));
                } else {
                     String hexString = colorSymbols.get(string.toLowerCase());
                     if(hexString != null) {
                         c = hexToColor(hexString);
                     }
                }


                if (c == null) {
                    c = Color.black;
                }
                colorCache.put(string, c);
            }
            return c;

        } catch (NumberFormatException numberFormatException) {
            log.error("Error in color string. ", numberFormatException);
            return Color.black;
        }
    }

    private static Color hexToColor(String string) {
        if (string.length() == 6) {
            int red = Integer.parseInt(string.substring(0, 2), 16);
            int green = Integer.parseInt(string.substring(2, 4), 16);
            int blue = Integer.parseInt(string.substring(4, 6), 16);
            return new Color(red, green, blue);
        } else {
            return null;
        }

    }



    /**
     * Return  alphas shaded color.  This method is used, rather than the Color constructor, so that
     * the alpha is not lost in postscript output.
     *
     * @param dest
     * @param source
     * @param alpha
     * @return
     */
    public static Color getCompositeColor(float[] dest, float[] source, float alpha) {
        int r = (int) ((alpha * source[0] + (1 - alpha) * dest[0]) * 255 + 0.5);
        int g = (int) ((alpha * source[1] + (1 - alpha) * dest[1]) * 255 + 0.5);
        int b = (int) ((alpha * source[2] + (1 - alpha) * dest[2]) * 255 + 0.5);
        int a = 255;
        int value = ((a & 0xFF) << 24) |
                ((r & 0xFF) << 16) |
                ((g & 0xFF) << 8) |
                ((b & 0xFF) << 0);

        Color c = colorCache.get(value);
        if (c == null) {
            c = new Color(value);
            colorCache.put(value, c);
        }
        return c;
    }

    /**
     * Return  alphas shaded color for a white background.  This method is used, rather than the Color constructor, so that
     * the alpha is not lost in postscript output.
     *
     * @param source
     * @param alpha
     * @return
     */
    public static Color getCompositeColor(float[] source, float alpha) {
        return getCompositeColor(whiteComponents, source, alpha);
    }
}
