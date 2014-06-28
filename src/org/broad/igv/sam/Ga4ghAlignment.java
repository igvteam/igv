package org.broad.igv.sam;

import org.broad.igv.feature.Strand;
import org.broad.igv.track.WindowFunction;

import java.awt.*;

/**
 * Created by jrobinso on 6/17/14.
 */
public class Ga4ghAlignment implements Alignment {
    @Override
    public String getReadName() {
        return null;
    }

    @Override
    public String getReadSequence() {
        return null;
    }

    @Override
    public String getChromosome() {
        return null;
    }

    @Override
    public String getChr() {
        return null;
    }

    @Override
    public int getStart() {
        return 0;
    }

    @Override
    public int getEnd() {
        return 0;
    }

    @Override
    public int getAlignmentStart() {
        return 0;
    }

    @Override
    public boolean contains(double location) {
        return false;
    }

    @Override
    public AlignmentBlock[] getAlignmentBlocks() {
        return new AlignmentBlock[0];
    }

    @Override
    public AlignmentBlock[] getInsertions() {
        return new AlignmentBlock[0];
    }

    @Override
    public char[] getGapTypes() {
        return new char[0];
    }

    @Override
    public String getCigarString() {
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
    public boolean isProperPair() {
        return false;
    }

    @Override
    public boolean isMapped() {
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
    public Strand getFirstOfPairStrand() {
        return null;
    }

    @Override
    public Strand getSecondOfPairStrand() {
        return null;
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
    public int getAlignmentEnd() {
        return 0;
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
    public String getSample() {
        return null;
    }

    @Override
    public String getReadGroup() {
        return null;
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
    public boolean isSmallInsert() {
        return false;
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
    public String getLibrary() {
        return null;
    }

    @Override
    public String getClipboardString(double location) {
        return null;
    }

    @Override
    public Strand getReadStrand() {
        return null;
    }

    @Override
    public void finish() {

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
    public String getValueString(double position, WindowFunction windowFunction) {
        return null;
    }
}
