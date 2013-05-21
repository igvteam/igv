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

import org.broad.igv.data.Interval;

/**
 * Represents a cufflinks value from any of a fpkm tracking file as described here
 *    http://cufflinks.cbcb.umd.edu/manual.html#fpkm_tracking_format
 * @author jrobinso
 *         Date: 3/8/13
 *         Time: 2:37 PM
 */
public class FPKMValue extends Interval {

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
        return new FPKMSampleValue(chr, start, end, chr, fpkm[sampleIndex], fpkmLo[sampleIndex], fpkmHi[sampleIndex]);
    }

    public String getGene(){
        return this.gene;
    }

    public int getNumSamples(){
        return this.fpkm.length;
    }



}
