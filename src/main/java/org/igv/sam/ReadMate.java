package org.igv.sam;


import htsjdk.samtools.util.Locatable;
import org.igv.feature.Strand;

public class ReadMate implements Locatable {

    private String chr;
    int start;
    private boolean negativeStrand;
    boolean mapped;

    public ReadMate(String chr, int start, boolean negativeStrand,
                    boolean isReadUnmappedFlag) {
        this.chr = chr;
        this.start = start;
        this.negativeStrand = negativeStrand;
        this.mapped = !isReadUnmappedFlag && !chr.equals("*");
    }

    public boolean isMapped() {
        return mapped;
    }

    public String positionString() {
        return chr + ":" + start + " (" + (isNegativeStrand() ? "-" : "+") + ")";
    }

    @Override
    public String getContig() {
        return chr;
    }

    public int getStart() {
        return start;
    }

    @Override
    public int getEnd() {
        throw new UnsupportedOperationException("ReadMate has an unknown end.");
    }

    public boolean isNegativeStrand() {
        return negativeStrand;
    }

    public Strand getStrand() {
        return negativeStrand ? Strand.NEGATIVE : Strand.POSITIVE;
    }

    public String getChr() {
        return getContig();
    }
}
