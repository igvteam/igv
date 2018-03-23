/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
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
