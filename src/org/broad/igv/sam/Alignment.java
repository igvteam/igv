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

import java.awt.*;

/**
 * @author jrobinso
 */
public interface Alignment extends LocusScore {

    String getReadName();

    String getReadSequence();

    String getChr();

    int getAlignmentStart();

    int getAlignmentEnd();

    boolean contains(double location);

    AlignmentBlock[] getAlignmentBlocks();

    AlignmentBlock[] getInsertions();

    String getCigarString();

    char[] getGapTypes();

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
    Color getColor();

    String getSample();

    String getReadGroup();

    String getLibrary();

    String getClipboardString(double location);


    void finish();
}
