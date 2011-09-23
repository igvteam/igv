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
import java.util.Hashtable;
import java.util.Map;

/**
 * @author Jim Robinson
 * @date 9/22/11
 */
public class ColorScale {

    Map<Integer, Color> colorCache = new Hashtable();

    int minCount = 0;
    int maxCount = 200;
    Color background;

    public Color getColor(double score) {

        if(score < minCount) {
            return background;
        }

        float alpha = (float) Math.max(0.05f, Math.min(1.0f, score / maxCount));


        float[] comps = background.getColorComponents(null);

        int idx = (int) (100 * alpha);
        Color c = colorCache.get(idx);
        if (c == null) {
            float rAlpha = Math.max(0.05f, Math.min(1.0f, 0.01f * idx));
            float red = ((1 - rAlpha) * comps[0] + rAlpha);
            float green = ((1 - rAlpha) * comps[1]);
            float blue = ((1 - rAlpha) * comps[2]);
            c = new Color(red, green, blue);
            colorCache.put(idx, c);
        }
        return c;
    }
}
