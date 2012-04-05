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

/**
 * @author Neva Cherniavsky
 * @date 3/22/12
 */
public class HiCColorScale implements org.broad.igv.renderer.ColorScale {

    private float min = -1f;
    private float max = 1f;

    public HiCColorScale() {
    }

    public void setMin(float min) {
        this.min = min;
        System.out.println(min);
    }

    public void setMax(float max) {
        this.max = max;
        System.out.println(max);
    }

    public Color getColor(float score) {

        if(score > 0) {
            score = score/max;
            int R = (int) ( 255 * Math.min(score,1));
            int G = 0;
            int B = 0;
            return new Color(R,G,B);
        } else if(score < 0) {
            score = score/min;
            int R = 0;
            int G = 0;
            int B = (int) (255 * Math.min(score,1));
            return new Color(R,G,B);
        } else {
            // Nan ?
            return Color.black;
        }

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

}