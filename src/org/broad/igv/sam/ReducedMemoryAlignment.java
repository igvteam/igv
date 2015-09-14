/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.Strand;
import org.broad.igv.track.WindowFunction;

import java.awt.*;

/**
 * @author jrobinso
 */
public class ReducedMemoryAlignment implements Alignment {

    //String readName;
    private String chromosome;
    private int start;
    private int end;
    boolean negativeStrand;


    public ReducedMemoryAlignment(Alignment al) {
        this.negativeStrand = al.isNegativeStrand();
      //  this.readName = al.getReadName();
        this.chromosome = al.getChr();
        this.start = al.getStart();
        this.end = al.getEnd();
    }



    public String getReadName() {
        return null;
    }

    /**
     * .aligned files do not include sequence
     *
     * @return
     */
    public String getReadSequence() {
        return "";
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


    public Color getColor() {
        return null;
    }


    public String getChromosome() {
        return chromosome;
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

    public String getClipboardString(double location) {
        return getValueString(location, null);
    }

    public String getValueString(double position, WindowFunction windowFunction) {
        return "Read length = " + (getEnd() - getStart());
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
        return isNegativeStrand() ? Strand.NEGATIVE : Strand.POSITIVE;
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
