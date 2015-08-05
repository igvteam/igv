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

/*
 * AlphaColorGradient.java
 *
 * Created on November 17, 2007, 3:16 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package org.broad.igv.util;

import java.awt.*;

/**
 * @author jrobinso
 */
public class AlphaColorGradient {
    Color maxColor = Color.red;
    Color minColor = Color.blue;
    double minValue = -1.5;
    double midValue = 0;
    double maxValue = 1.5;

    Color getColor(double value) {
        if (value > midValue) {
            int alpha = (int) (255 * (value - midValue) / (maxValue - midValue));
            alpha = Math.max(0, Math.min(255, alpha));
            return new Color(maxColor.getRed(), maxColor.getGreen(), maxColor.getBlue(), alpha);
        } else {
            int alpha = (int) (255 * (midValue - value) / (midValue - minValue));
            alpha = Math.max(0, Math.min(255, alpha));
            return new Color(minColor.getRed(), minColor.getGreen(), minColor.getBlue(), alpha);
        }
    }
}

