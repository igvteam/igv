/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.igv.sam;

import org.igv.feature.LocusScore;
import org.igv.feature.Strand;
import org.igv.track.WindowFunction;

import java.awt.*;
import java.util.List;

/**
 * @author jrobinso
 */
public class DotAlignedAlignment implements Alignment {

    String readName;
    private String chromosome;
    private int start;
    private int end;
    boolean negativeStrand;


    public DotAlignedAlignment(String chromosome, int start, int end, boolean isNegative, String name) {
        this.negativeStrand = isNegative;
        this.readName = name;
        this.chromosome = chromosome;
        this.start = start;
        this.end = end;
    }

    public DotAlignedAlignment(String chromosome, int start, int end, boolean negativeStrand) {
        this.readName = chromosome + ":" + (start + 1) + "-" + (end + 1) + "(" +
                (negativeStrand ? "-" : "+") + ")";
        this.chromosome = chromosome;
        this.start = start;
        this.end = end;
        this.negativeStrand = negativeStrand;
    }

    public String getReadName() {
        return readName;
    }

    public void setMateSequence(String sequnce) {
        // Ignore
    }

    public String getPairOrientation() {
        return "";
    }

    public boolean isSmallInsert() {
        return false;
    }

    public boolean isVendorFailedRead() {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }


    public Color getYcColor() {
        return null;
    }

    public String getChr() {
        return chromosome;
    }

    @Override
    public String getContig() {
        return chromosome;
    }

    public int getAlignmentStart() {
        return getStart();
    }

    public boolean contains(double location) {
        return location >= getStart() && location < getEnd();
    }

    public AlignmentBlock[] getAlignmentBlocks() {
        return null;
    }

    public AlignmentBlock[] getInsertions() {
        return null;
    }

    public int getInferredInsertSize() {
        return 0;
    }

    public int getMappingQuality() {
        return 255;
    }

    public boolean isProperPair() {
        return true;
    }

    public boolean isMapped() {
        return true;
    }

    public boolean isPaired() {
        return false;
    }

    public boolean isNegativeStrand() {
        return negativeStrand;
    }

    public boolean isDuplicate() {
        return false;
    }

    public float getScore() {
        return 1.0f;
    }

    public LocusScore copy() {
        return this;
    }

    public String getValueString(double position, int mouseX, WindowFunction windowFunction) {
        return readName + "<br>Read length = " + (getEnd() - getStart());
    }

    /**
     * @return the start
     */
    public int getStart() {
        return start;
    }

    /**
     * @param start the start to set
     */
    public void setStart(int start) {
        this.start = start;
    }

    /**
     * @return the end
     */
    public int getEnd() {
        return end;
    }

    public int getAlignmentEnd() {
        return end;
    }

    /**
     * @param end the end to set
     */
    public void setEnd(int end) {
        this.end = end;
    }

    public byte getBase(double position) {
        return 0;
    }

    public byte getPhred(double position) {
        return 0;
    }

    @Override
    public List<Gap> getGaps() {
        return null;
    }


    public boolean isFirstOfPair() {
        return false;
    }

    public boolean isSecondOfPair() {
        return false;
    }

    public Strand getFirstOfPairStrand() {
        return isNegativeStrand() ? Strand.NEGATIVE : Strand.POSITIVE;
    }


    public Strand getSecondOfPairStrand() {
        return Strand.NONE;
    }

    public Strand getReadStrand() {
        return isNegativeStrand() ? Strand.NEGATIVE : Strand.POSITIVE;
    }

    @Override
    public boolean isPrimary() {
        return true;
    }

    @Override
    public boolean isSupplementary() {
        return false;
    }
}
