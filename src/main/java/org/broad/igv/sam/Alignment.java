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
import org.broad.igv.sam.mods.BaseModificationUtils;
import org.broad.igv.sam.mods.BaseModificationSet;
import org.broad.igv.track.WindowFunction;

import java.awt.*;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 */
public interface Alignment extends LocusScore {

    /**
     * Name for the associated read.  Cannot be null.
     * @return
     */
    default String getReadName() {
        return "";
    }

    default String getReadSequence() {
        return "";
    }

    String getChr();

    int getAlignmentStart();

    int getAlignmentEnd();

    boolean contains(double location);

    AlignmentBlock[] getAlignmentBlocks();

    AlignmentBlock[] getInsertions();

    String getCigarString();

    java.util.List<Gap> getGaps();

    int getInferredInsertSize();

    int getMappingQuality();

    ReadMate getMate();

    Strand getReadStrand();

    boolean isProperPair();

    boolean isMapped();

    boolean isPaired();

    boolean isFirstOfPair(); // Ben Berman

    boolean isSecondOfPair(); // Ben Berman

    boolean isNegativeStrand();

    boolean isDuplicate();

    boolean isPrimary();

    boolean isSupplementary();

    byte getBase(double position);

    byte getPhred(double position);

    Object getAttribute(String key);

    void setMateSequence(String sequence);

    String getPairOrientation();

    Strand getFirstOfPairStrand();

    Strand getSecondOfPairStrand();

    boolean isVendorFailedRead();

    /**
     * Return an explicitly set color for this alignment, if any  (typically null).
     * @return
     */
    Color getYcColor();

    String getSample();

    String getReadGroup();

    String getLibrary();

    String getClipboardString(double location, int mouseX);

    void finish();

    default AlignmentBlock getInsertionAt(int position) {
        final AlignmentBlock[] insertions = getInsertions();
        if(insertions == null) {
            return null;
        }
        for (AlignmentBlock block : insertions) {
            if (block.getStart() == position) return block;
        }
        return null;
    }

     default Gap getDeletionAt(int position) {
         List<Gap> gaps = this.getGaps();
         if (gaps != null && !gaps.isEmpty()) {
             for (Gap gap : gaps) {
                 if (gap.getStart() <= position
                         && gap.getnBases() + gap.getStart() > position
                         && gap.getType() == SAMAlignment.DELETION) {
                     return gap;
                 }
             }
         }
         return null;
     }

    default void setHaplotypeName(String hap) {}

    default String getHaplotypeName() {return null;}

    default void setHapDistance(int dist) {};

    default int getHapDistance() {return 0;}

    default Map<Integer, BaseModificationUtils> getBaseModificationMap() { return null;}

    default List<BaseModificationSet> getBaseModificationSets() { return null;}

    default short[] getSmrtSubreadIpd() { return null;}

    default short[] getSmrtSubreadPw() { return null;}

    default short[] getSmrtCcsIpd(boolean isForwardStrand) { return null;}

    default short[] getSmrtCcsPw(boolean isForwardStrand) { return null;}

    default String getAlignmentValueString(double position, int mouseX, AlignmentTrack.RenderOptions renderOptions) {
        return getValueString(position, mouseX, (WindowFunction) null);
    }

    /**
     * Get the most specific sub alignment which contains the given chromosome location, implementations which
     * contain multiple distinct sub alignments should override this to provide the appropriate behavior
     * Note: no check is performed to validate that the location is on the same chromosome as this alignment
     * @param location the location on the chromosome to select an alignment from
     * @return the alignment (or sub-alignment) that contains the given location, null if this alignment does not contain it
     */
    default Alignment getSpecificAlignment(double location) {
        return this.contains(location) ? this : null;
    }
}
