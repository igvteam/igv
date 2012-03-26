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

package org.broad.igv.hic;

import java.awt.*;

/**
 * @author Neva Cherniavsky
 * @date 3/22/12
 */
public class HiCColorScale implements org.broad.igv.renderer.ColorScale {

    // Use enums to avoid mistakes in initialization
    enum Scheme {
        ZERO, ONE, MINUS_ONE
    }

    private int maxval = 5; // equivalent to maxdiv = 5 in Erez's Python code
    private int minval = 0;
    private Scheme scheme = Scheme.ONE;
    private Color noDataColor = Color.BLACK;

    // Create a cache to limit the number of distinct colors created.  This is particularly important for performance
    // on linux/unix machines running X-Windows.
    private int nColors = 100;
    private float step;
    Color[] cachedColors;

    public HiCColorScale(Scheme scheme, int minval, int maxval) {
        this.maxval = maxval;
        this.minval = minval;
        this.scheme = scheme;
    }


    private void initColors() {
        cachedColors = new Color[nColors];
        step = ((float) (maxval - minval)) / (nColors - 1);
        for (int i = 0; i < nColors; i++) {
            float score = minval + i*step;
            cachedColors[i] = computeColor(score);
        }

    }

    public Color getColor(float score) {
        if (cachedColors == null) {
            initColors();
        }
        int idx = (int) ((score - minval) / step);
        if (idx < 0) {
            return cachedColors[0];
        } else if (idx >= cachedColors.length) {
            return cachedColors[cachedColors.length - 1];
        } else return cachedColors[idx];
    }

    private Color computeColor(float score) {
        int red, green, blue;

        if(Float.isNaN(score)) {
            return noDataColor;
        }

        if (scheme == Scheme.ONE) {
            if (score <= 0) {
                return Color.BLACK;
            }
            red = (int) (255 * Math.min(score / maxval, 1));
            green = 0;
            blue = (int) (255 * Math.min(Math.pow(score, -1) / maxval, 1));
            return new Color(red, green, blue);

        } else if (scheme == Scheme.ZERO) {
            red = 255;
            green = (int) (Math.max(0, 255 - 255 * (score - minval) / (maxval - minval)));
            blue = (int) (Math.max(0, 255 - 255 * (score - minval) / (maxval - minval)));
            if (green < 0) green = 0;
            if (blue < 0) blue = 0;
            if (green > 255) green = 255;
            if (blue > 255) blue = 255;
            return new Color(red, green, blue);

        } else if (scheme == Scheme.MINUS_ONE) {
            if (score > 0) {
                red = (int) (255 * Math.min(score / (maxval), 1));
                green = 0;
                blue = 0;
                return new Color(red, green, blue);
            } else if (score < 0) {
                red = 0;
                green = 0;
                blue = (int) (255 * Math.min(((-1) * score) / (maxval), 1));
                return new Color(red, green, blue);
            } else {
                return Color.black;
            }
        } else return Color.black;
    }

    public Color getColor(String symbol) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }


    public Color getNoDataColor() {
        return noDataColor;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public String asString() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public boolean isDefault() {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

}