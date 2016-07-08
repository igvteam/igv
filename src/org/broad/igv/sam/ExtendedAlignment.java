package org.broad.igv.sam;

import org.broad.igv.feature.Strand;
import org.broad.igv.track.WindowFunction;

import java.awt.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Class for experimenting with barcodes
 */
public class ExtendedAlignment implements Alignment {


    final String barcode;
    List<Alignment> alignments;
    int alignmentStart;
    int alignmentEnd;
    String chr;

    public ExtendedAlignment(String barcode) {
        this.barcode = barcode;
        alignments = new ArrayList<>();
    }

    public void addAlignment(Alignment alignment) {

        if (alignments.isEmpty()) {
            this.chr = alignment.getChr();
            this.alignmentStart = alignment.getAlignmentStart();
            this.alignmentEnd = alignment.getAlignmentEnd();
        } else {
            if (!this.chr.equals(alignment.getChr())) {
                throw new RuntimeException("Mixed chromosome extended alignments not supported");
            }
            this.alignmentStart = Math.min(alignment.getAlignmentStart(), this.alignmentStart);
            this.alignmentEnd = Math.max(alignment.getAlignmentEnd(), this.alignmentEnd);
        }

        this.alignments.add(alignment);
    }

    public String getBarcode() {
        return barcode;
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

        buffer.append("Barcode: " + this.barcode);

        for(Alignment a : alignments) {
            if(a.contains(position)) {
                buffer.append("----------------------<br>");
                buffer.append(a.getValueString(position, windowFunction));
                break;
            }
        }

        return buffer.toString();
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
    public int getMappingQuality() {
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
    public Object getAttribute(String key) {
        return null;
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
