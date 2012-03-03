/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
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


    public Color getDefaultColor() {
        return AlignmentRenderer.grey1;
    }


    public String getChromosome() {
        return chromosome;
    }

    public String getChr() {
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
}
