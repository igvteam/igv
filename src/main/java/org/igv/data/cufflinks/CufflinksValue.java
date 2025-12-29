package org.igv.data.cufflinks;

import org.igv.feature.LocusScore;
import org.igv.feature.Range;

/**
 * @author jrobinso
 *         Date: 3/8/13
 *         Time: 9:32 PM
 */
abstract public class CufflinksValue extends Range implements LocusScore {
    String gene;

    public CufflinksValue(String chr, int start, int end, String gene) {
        super(chr, start, end);
        this.gene = gene;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public String getGene() {
        return gene;
    }

}
