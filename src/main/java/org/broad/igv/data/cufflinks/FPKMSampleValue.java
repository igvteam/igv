package org.broad.igv.data.cufflinks;

import org.broad.igv.track.WindowFunction;

/**
 * Represents a value from a cufflinks file for a single sample
 * @see FPKMValue
 * @author jrobinso, jacob
 */
public class FPKMSampleValue extends CufflinksValue{

    float fpkm;
    float fpkmLo;
    float fpkmHi;

    public FPKMSampleValue(String chr, int start, int end, String gene, float fpkm, float fpkmLo, float fpkmHi) {
        super(chr, start, end, gene);
        this.fpkm = fpkm;
        this.fpkmLo = fpkmLo;
        this.fpkmHi = fpkmHi;
    }

    @Override
    public float getScore() {
        return fpkm;
    }

    @Override
    public String getValueString(double position, int mouseX, WindowFunction windowFunction) {

        StringBuilder sb = new StringBuilder();
        sb.append(getChr() + ":" + (getStart() + 1) + "-" + getEnd());
        sb.append("<br>Gene = " + gene);
        sb.append("<br>FPKM = " + fpkm);
        sb.append("<br>FPKM_LO = " + fpkmLo);
        sb.append("<br>FPKM_HI = " + fpkmHi);
        return sb.toString();
    }

}
