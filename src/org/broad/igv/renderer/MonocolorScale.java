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

package org.broad.igv.renderer;

import java.awt.*;

/**
 * @author jrobinso
 * @date Jul 11, 2011
 */
public class MonocolorScale extends AbstractColorScale {

    static int nSteps = 10;
    float maxAlpha = 1.0f;
    float minAlpha = 0.1f;

    float minValue;
    float maxValue;
    Color [] colorTable = new Color[nSteps];
    private final float range;


    public MonocolorScale(float x1, float x2, Color baseColor) {

        minValue = Math.min(x1, x2);
        maxValue = Math.max(x1, x2);
        range = maxValue - minValue;

        if(x1 == x2) {
            for(int i=0; i<colorTable.length; i++)  {
                colorTable[i] = baseColor;
            }
        }
        else {
            float alphaStep = (maxAlpha - minAlpha) / (nSteps - 1);
            float [] baseComponents = baseColor.getComponents(null);
            for(int i=0; i<nSteps; i++) {
                float alpha =  Math.min(1, 0.1f + minAlpha + i*alphaStep);
                colorTable[i] = new Color(baseComponents[0], baseComponents[1], baseComponents[2], alpha);
            }
        }

    }




    @Override
    public Color getColor(float value) {
        if(range == 0) {
            return colorTable[0];
        }

        int idx = (int) Math.min(nSteps - 1, Math.max(0, nSteps * (value - minValue) / range));
        return colorTable[idx];
    }


    // TODO -- implement this so the class can be serialized
    public String asString() {
        return null;
    }

    public boolean isDefault() {
        return false;  
    }
}
