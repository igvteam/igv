package org.broad.igv.sam;

import org.broad.igv.Globals;
import org.broad.igv.feature.Strand;
import org.broad.igv.track.WindowFunction;

import java.awt.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Class for experimenting with 10X linked reads.
 */

public class LinkedAlignment implements Alignment {


    final String tag;
    final String barcode;
    String haplotype;
    Strand strand;
    List<Alignment> alignments;
    String chr;
    int alignmentStart;
    int alignmentEnd;
    int lastAlignmentEnd = 0;
    int maxGap = 0;
    Map<String, Object> attributes;

    public LinkedAlignment(String tag, String bc) {
        attributes = new HashMap<>();
        this.tag = tag;
        this.barcode = bc;
        attributes.put(tag, barcode);
        alignments = new ArrayList<>();
    }

    public void addAlignment(Alignment alignment) {

        if (alignments.isEmpty()) {
            this.chr = alignment.getChr();
            alignmentStart = alignment.getAlignmentStart();
            alignmentEnd = alignment.getAlignmentEnd();
            lastAlignmentEnd = alignment.getAlignmentEnd();

            Object hp = alignment.getAttribute("HP");
            haplotype = hp == null ? null : hp.toString();

            strand = alignment.getReadStrand();

        } else {

            if (!this.chr.equals(alignment.getChr())) {
                throw new RuntimeException("Mixed chromosome linked alignments not supported");
            }
            maxGap = Math.max(alignment.getAlignmentStart() - lastAlignmentEnd, maxGap);
            lastAlignmentEnd = alignment.getAlignmentEnd();
            alignmentEnd = Math.max(alignment.getAlignmentEnd(), this.alignmentEnd);

            Object hp = alignment.getAttribute("HP");
            if (hp != null) {
                if (!hp.toString().equals(this.haplotype)) {
                    this.haplotype = "MIXED";
                }
            }


            if (strand != alignment.getReadStrand()) {
                strand = Strand.NONE;   // i.e. mixed
            }
        }

        alignments.add(alignment);
    }

    public String getBarcode() {
        return barcode;
    }

    public Strand getStrand() {
        return strand;
    }

    @Override
    public String getChr() {
        return chr;
    }

    @Override
    public int getAlignmentStart() {
        return this.alignmentStart;
    }

    @Override
    public int getAlignmentEnd() {
        return this.alignmentEnd;
    }


    @Override
    public int getStart() {
        return this.alignmentStart;

    }

    @Override
    public int getEnd() {
        return this.alignmentEnd;
    }

    @Override
    public boolean contains(double location) {
        return location >= this.alignmentStart && location <= this.alignmentEnd;
    }

    @Override
    public boolean isMapped() {
        return true;
    }

    @Override
    public String getValueString(double position, WindowFunction windowFunction) {
        StringBuffer buffer = new StringBuffer();

        buffer.append("Linking id (" + tag + ") = " + this.barcode);
        if (this.haplotype != null) buffer.append("<br>Haplotype = " + this.haplotype);
        buffer.append("<br># alignments = " + alignments.size());
        buffer.append("<br>Total span = " + Globals.DECIMAL_FORMAT.format(getAlignmentEnd() - getAlignmentStart()) + "bp<br>");


        for (Alignment a : alignments) {
            if (a.contains(position)) {
                buffer.append("----------------------<br>");
                buffer.append(a.getValueString(position, windowFunction));
                break;
            }
        }

        return buffer.toString();
    }


    @Override
    public Object getAttribute(String key) {
        if ("HP".equals(key)) {
            return haplotype;
        } else {
            return attributes.get(key);
        }
    }


    @Override
    public int getMappingQuality() {
        return 30;    // This is used for coloring.  Not sure what to do here
    }

    /////////////////////////////////////////////////////////////

    @Override
    public String getReadName() {
        return null;
    }

    @Override
    public String getReadSequence() {
        return null;
    }

    @Override
    public AlignmentBlock[] getAlignmentBlocks() {
        return null;
    }

    @Override
    public AlignmentBlock[] getInsertions() {
        return null;
    }

    @Override
    public String getCigarString() {
        return null;
    }

    @Override
    public List<Gap> getGaps() {
        return null;
    }

    @Override
    public int getInferredInsertSize() {
        return 0;
    }

    @Override
    public ReadMate getMate() {
        return null;
    }

    @Override
    public Strand getReadStrand() {
        return null;
    }

    @Override
    public boolean isProperPair() {
        return false;
    }


    @Override
    public boolean isPaired() {
        return false;
    }

    @Override
    public boolean isFirstOfPair() {
        return false;
    }

    @Override
    public boolean isSecondOfPair() {
        return false;
    }

    @Override
    public boolean isNegativeStrand() {
        return false;
    }

    @Override
    public boolean isDuplicate() {
        return false;
    }

    @Override
    public boolean isPrimary() {
        return false;
    }

    @Override
    public boolean isSupplementary() {
        return false;
    }

    @Override
    public byte getBase(double position) {
        return 0;
    }

    @Override
    public byte getPhred(double position) {
        return 0;
    }

    @Override
    public void setMateSequence(String sequence) {

    }

    @Override
    public String getPairOrientation() {
        return null;
    }

    @Override
    public Strand getFirstOfPairStrand() {
        return null;
    }

    @Override
    public Strand getSecondOfPairStrand() {
        return null;
    }

    @Override
    public boolean isVendorFailedRead() {
        return false;
    }

    @Override
    public Color getColor() {
        return null;
    }

    @Override
    public String getSample() {
        return null;
    }

    @Override
    public String getReadGroup() {
        return null;
    }

    @Override
    public String getLibrary() {
        return null;
    }

    @Override
    public String getClipboardString(double location) {
        return null;
    }

    @Override
    public void finish() {

    }

    @Override
    public void setStart(int start) {

    }

    @Override
    public void setEnd(int end) {

    }

    @Override
    public float getScore() {
        return 0;
    }

    @Override
    public String getContig() {
        return null;
    }

}
