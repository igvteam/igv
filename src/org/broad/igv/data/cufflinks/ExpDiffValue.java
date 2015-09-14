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

package org.broad.igv.data.cufflinks;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.WindowFunction;

/**
 * @author jrobinso
 *         Date: 3/8/13
 *         Time: 9:33 PM
 */
public class ExpDiffValue extends CufflinksValue implements LocusScore {


    float log2Ratio;
    float fpkmX;
    float fpkmY;
    String significant;

    public ExpDiffValue(String chr, int start, int end, String gene, float log2Ratio, float fpkmX, float fpkmY, String significant) {
        super(chr, start, end, gene);
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
