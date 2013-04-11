/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
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
 * @author jrobinso
 *         Date: 3/8/13
 *         Time: 9:33 PM
 */
public class ExpDiffValue extends CufflinksValue {


    float log2Ratio;
    float fpkmX;
    float fpkmY;
    String significant;

    public ExpDiffValue(String gene, String chr, int start, int end, float log2Ratio, float fpkmX, float fpkmY, String significant) {
        super(gene, chr, start, end);
        this.log2Ratio = log2Ratio;
        this.fpkmX = fpkmX;
        this.fpkmY = fpkmY;
        this.significant = significant;
    }

    @Override
    public float getScore() {
       return log2Ratio;
    }

    @Override
    public String getValueString(double position, WindowFunction windowFunction) {

        StringBuilder sb = new StringBuilder();
        sb.append(getChr() + ":" + (getStart() + 1) + "-" + getEnd());
        sb.append("<br>gene = " + gene);
        sb.append("<br>log2(y/x) = " + log2Ratio);
        sb.append("<br>FPKM X = " + fpkmX);
        sb.append("<br>FPKM Y = " + fpkmY);
        sb.append("<br>Significant? " + significant);
        return sb.toString();   }
}
