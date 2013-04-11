/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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

    float fpkm;
    float fpkmLo;
    float fpkmHi;

    public FPKMValue(String gene, String chr, int start, int end, float fpkm, float fpkmLo, float fpkmHi) {
        super(gene, chr, start, end);
        this.fpkm = fpkm;
        this.fpkmLo = fpkmLo;
        this.fpkmHi = fpkmHi;
    }

    @Override
    public float getScore() {
        return fpkm;
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

    public String getGene() {
        return gene;
    }
}
