package org.igv.data.cufflinks;

import org.igv.feature.LocusScore;
import org.igv.track.WindowFunction;

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
    public String getValueString(double position, int mouseX, WindowFunction windowFunction) {

        StringBuilder sb = new StringBuilder();
        sb.append(getChr() + ":" + (getStart() + 1) + "-" + getEnd());
        sb.append("<br>gene = " + gene);
        sb.append("<br>log2(y/x) = " + log2Ratio);
        sb.append("<br>FPKM X = " + fpkmX);
        sb.append("<br>FPKM Y = " + fpkmY);
        sb.append("<br>Significant? " + significant);
        return sb.toString();   }
}
