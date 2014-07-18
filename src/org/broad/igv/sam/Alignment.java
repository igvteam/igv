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

    boolean isSmallInsert();

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
