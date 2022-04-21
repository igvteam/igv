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
import java.util.Map;

/**
 * @author jrobinso
 */
public interface Alignment extends LocusScore {

    String getReadName();

    String getChr();

    int getAlignmentStart();

    int getAlignmentEnd();

    Strand getReadStrand();

    String getClipboardString(double location, int mouseX);

    default boolean isNegativeStrand() {
        return Strand.NEGATIVE == getReadStrand();
    }

    default String getReadSequence() {
        return "*";
    }

    default int getReadSequenceLength() {
        return getReadSequence().length();
    }

    default boolean contains(double location) {
        return location >= getStart() && location < getEnd();
    }

    default AlignmentBlock[] getAlignmentBlocks() {
        return null;
    }

    default AlignmentBlock[] getInsertions() {
        return null;
    }

    default String getCigarString() {
        return "*";
    }

    default java.util.List<Gap> getGaps() {
        return null;
    }

    default int getInferredInsertSize() {
        return 0;
    }

    default int getMappingQuality() {
        return 255;
    }

    default ReadMate getMate() {
        return null;
    }

    default boolean isProperPair() {
        return true;
    }

    default boolean isMapped() {
        return true;
    }

    default boolean isPaired() {
        return false;
    }

    default boolean isFirstOfPair() {
        return false;
    }

    default boolean isSecondOfPair() {
        return false;
    }

    default boolean isDuplicate() {
        return false;
    }

    default boolean isPrimary() {
        return true;
    }

    default boolean isSupplementary() {
        return false;
    }

    default byte getBase(double position) {
        return 0;
    }

    default byte getPhred(double position) {
        return 0;
    }

    default Object getAttribute(String key) {
        return null;
    }

    default void setMateSequence(String sequence) {}

    default String getPairOrientation() {
        return null;
    }

    default Strand getFirstOfPairStrand() {
        return Strand.NONE;
    }

    default Strand getSecondOfPairStrand() {
        return Strand.NONE;
    }

    default boolean isVendorFailedRead() {
        return false;
    }

    /**
     * Return an explicitly set color for this alignment, if any  (typically null).
     *
     * @return
     */
    default Color getYcColor() {
        return null;
    }

    default String getSample() {
        return null;
    }

    default String getReadGroup() {
        return null;
    }

    default String getLibrary() {
        return null;
    }

    default AlignmentBlock getInsertionAt(int position) {
        return null;
    }

    default Map<Integer, BaseModification> getBaseModificationMap() {
        return null;
    }

    default String getAlignmentValueString(double position, int mouseX, AlignmentTrack.RenderOptions renderOptions) {
        return getValueString(position, mouseX, (WindowFunction) null);
    }

    default void finish() {}


    // Experimental
    default void setHaplotypeName(String hap) {
    }

    default String getHaplotypeName() {
        return null;
    }

}
