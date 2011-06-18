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

