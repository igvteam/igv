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

import org.broad.igv.renderer.ContinuousColorScale;

import java.awt.*;
import java.util.Hashtable;
import java.util.Map;

/**
 * @author Jim Robinson
 * @date 9/22/11
 */
public class ObservedColorScale implements org.broad.igv.renderer.ColorScale {

    private int colorTableSize = 100;
    Color [] colorCache = new Color[colorTableSize];  // Array instead of map for efficiency

    private int minCount = 0;
    private int maxCount = 200;
    private Color background;
    private float[] backgroundColorComponents;


    public Color getColor(float score) {

        if (score < minCount) {
            return background;
        }

        final float v = score / maxCount;
        float alpha = v < 0.05f ? 0.05f : (v > 1.0f ? 1.0f : v);   // Math.min and Math.max taking significant time
        int idx = (int) ((colorTableSize - 1) * alpha);
        Color c = colorCache[idx];
        if (c == null) {
            float rAlpha = Math.max(0.05f, Math.min(1.0f, 0.01f * idx));
            float red = ((1 - rAlpha) * backgroundColorComponents[0] + rAlpha);
            float green = ((1 - rAlpha) * backgroundColorComponents[1]);
            float blue = ((1 - rAlpha) * backgroundColorComponents[2]);
            c = new Color(red, green, blue);
            colorCache[idx] = c;
        }
        return c;
    }

    public Color getColor(String symbol) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }


    public Color getNoDataColor() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public String asString() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public boolean isDefault() {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public int getMinCount() {
        return minCount;
    }

    public int getMaxCount() {
        return maxCount;
    }

    public Color getBackground() {
        return background;
    }

    public void setMinCount(int minCount) {
        this.minCount = minCount;
    }

    public void setMaxCount(int maxCount) {
        this.maxCount = maxCount;
    }

    public void setBackground(Color background) {
        this.background = background;
        backgroundColorComponents = background.getColorComponents(null);
    }

    public void setRange(int min, int max) {
        this.minCount = min;
        this.maxCount = max;

    }
}
