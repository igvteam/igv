package org.igv.data.cufflinks;

import org.igv.feature.Range;

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
