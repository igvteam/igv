/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.sam;

import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.Exon;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.Strand;
import org.broad.igv.track.WindowFunction;

import java.awt.*;

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
    AlignmentBlock[] blocks;
    Strand strand;

    public FeatureWrappedAlignment(BasicFeature f) {

        this.readName = f.getName();
        this.chr = f.getChr();
        this.start = f.getStart();
        this.end = f.getEnd();
        strand = f.getStrand();

        if (f.getExonCount() > 0) {
            blocks = new AlignmentBlock[f.getExonCount()];
            int i = 0;
            for (Exon exon : f.getExons()) {
                int length = exon.getLength();
                byte[] seq = new byte[length];
                blocks[i] = new AlignmentBlock(getChr(), exon.getStart(), seq, seq);
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

    public String getCigarString() {
        return "*";
    }

    public int getInferredInsertSize() {
        return 0;
    }

    public int getMappingQuality() {
        return 255;
    }

    public ReadMate getMate() {
        return null;
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

    public String getClipboardString(double location) {
        return getValueString(location, null);
    }

    public String getValueString(double position, WindowFunction windowFunction) {
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

    public String getSample() {
        return null;
    }

    public String getReadGroup() {
        return null;
    }

    public String getLibrary() {
        return null;
    }

    public Object getAttribute(String key) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
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

    public Color getColor() {
        return null;
    }

    public char[] getGapTypes() {
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
    public void finish() {
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
