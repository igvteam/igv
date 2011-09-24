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

import javax.swing.*;
import java.awt.*;
import java.io.Serializable;

/**
 * @author Jim Robinson
 * @date 9/22/11
 */
public class ColorScalePanel extends JComponent implements Serializable {

    ColorScale colorScale;

    public ColorScalePanel() {

    }

    public ColorScalePanel(ColorScale colorScale) {
        this.colorScale = colorScale;
    }


    @Override
    protected void paintComponent(Graphics graphics) {

        if (colorScale != null) {

            int nSteps = getWidth() - 1;
            if(nSteps <= 0) {
                return;
            }
            float delta = (colorScale.maxCount / nSteps);

            int xLast = 0;
            for (int n = 1; n < nSteps; n++) {
                int x = n;
                int dx = x - xLast;
                xLast = x;

                float v = delta * n;
                Color c = colorScale.getColor(v);
                graphics.setColor(c);
                graphics.fillRect(x, 0, dx, getHeight());

            }


        }
    }
}
