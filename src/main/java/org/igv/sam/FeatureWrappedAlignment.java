package org.igv.sam;

import org.igv.feature.BasicFeature;
import org.igv.feature.Exon;
import org.igv.feature.LocusScore;
import org.igv.feature.Strand;
import org.igv.track.WindowFunction;

import java.awt.*;
import java.util.List;

/**
 * Some alignment formats are parsed as Features.
 * <p/>
 * This is all getting rather circular,  some refactoring is in order.
 *
 * @author jrobinso
 * @date Aug 5, 2010
 */
public class FeatureWrappedAlignment implements Alignment {
    String readName;
    String chr;
    int start;
    int end;
    AlignmentBlockImpl[] blocks;
    Strand strand;

    public FeatureWrappedAlignment(BasicFeature f) {

        this.readName = f.getName() == null ? "" : f.getName();
        this.chr = f.getChr();
        this.start = f.getStart();
        this.end = f.getEnd();
        strand = f.getStrand();

        if (f.getExonCount() > 0) {
            blocks = new AlignmentBlockImpl[f.getExonCount()];
            int i = 0;
            for (Exon exon : f.getExons()) {
                int length = exon.getLength();
                byte[] seq = new byte[length];
                blocks[i] = new AlignmentBlockImpl(exon.getStart(), seq, seq);
                i++;
            }
        }
    }

    public String getReadName() {
        return readName;
    }

    public String getReadSequence() {
        return null;
    }

    public String getChromosome() {
        return chr;
    }

    public String getChr() {
        return chr;
    }

    @Override
    public String getContig() {
        return chr;
    }

    public int getAlignmentStart() {
        return start;
    }

    public boolean contains(double location) {
        return location >= start && location <= getEnd();
    }

    public AlignmentBlock[] getAlignmentBlocks() {
        return blocks;
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
        return strand == Strand.NEGATIVE;
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

    public void setMateSequence(String sequence) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public String getPairOrientation() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public boolean isSmallInsert() {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public boolean isVendorFailedRead() {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public Color getYcColor() {
        return null;
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
        return strand;
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
