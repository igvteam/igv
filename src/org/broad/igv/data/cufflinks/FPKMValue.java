/*
 * Copyright (c) 2007-2013 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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

package org.broad.igv.data.cufflinks;

import org.broad.igv.track.WindowFunction;

/**
 * Represents a cufflinks value from any of a fpkm tracking file as described here
 *    http://cufflinks.cbcb.umd.edu/manual.html#fpkm_track
 * @author jrobinso
 *         Date: 3/8/13
 *         Time: 2:37 PM
 */
public class FPKMValue extends CufflinksValue {

    String gene;
    float fpkm;
    float fpkmLo;
    float fpkmHi;

    public FPKMValue(String gene, String chr, int start, int end, float fpkm, float fpkmLo, float fpkmHi) {
        super(chr, start, end);
        this.gene = gene;
        this.fpkm = fpkm;
        this.fpkmLo = fpkmLo;
        this.fpkmHi = fpkmHi;
    }

    @Override
    public float getScore() {
        return fpkm;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getValueString(double position, WindowFunction windowFunction) {

        StringBuilder sb = new StringBuilder();
        sb.append(getChr() + ":" + (getStart() + 1) + "-" + getEnd());
        sb.append("<br>Gene = " + gene);
        sb.append("<br>FPKM = " + fpkm);
        sb.append("<br>FPKM_LO = " + fpkmLo);
        sb.append("<br>FPKM_HI = " + fpkmHi);
        return sb.toString();
    }
}
