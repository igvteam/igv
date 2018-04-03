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

import org.broad.igv.feature.Range;

/**
 * Represents a cufflinks value from any of a fpkm tracking file as described here
 *    http://cufflinks.cbcb.umd.edu/manual.html#fpkm_tracking_format
 * @author jrobinso
 *         Date: 3/8/13
 *         Time: 2:37 PM
 */
public class FPKMValue extends Range {

    float[] fpkm;
    float[] fpkmLo;
    float[] fpkmHi;

    String gene;

    public FPKMValue(String chr, int start, int end, String gene, float[] fpkm, float[] fpkmLo, float[] fpkmHi) {
        super(chr, start, end);
        this.gene = gene;
        assert fpkm.length == fpkmLo.length && fpkm.length == fpkmHi.length;
        this.fpkm = fpkm;
        this.fpkmLo = fpkmLo;
        this.fpkmHi = fpkmHi;
    }

    public FPKMSampleValue getSampleValue(int sampleIndex){
        return new FPKMSampleValue(chr, start, end, gene, fpkm[sampleIndex], fpkmLo[sampleIndex], fpkmHi[sampleIndex]);
    }

    public String getGene(){
        return this.gene;
    }

    public int getNumSamples(){
        return this.fpkm.length;
    }



}
