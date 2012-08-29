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

import org.broad.tribble.Feature;

/**
 * @author Jim Robinson
 * @date 11/22/11
 */
public interface AlignmentCounts extends Feature {

    int getTotalCount(int pos);

    int getNegTotal(int pos);

    int getPosTotal(int pos);

    int getTotalQuality(int pos);

    int getCount(int pos, byte b);

    int getNegCount(int pos, byte b);

    int getPosCount(int pos, byte b);

    int getDelCount(int pos);

    int getInsCount(int pos);

    int getQuality(int pos, byte b);

    int getAvgQuality(int pos, byte b);

    void incCounts(Alignment alignment);

    int getNumberOfPoints();

    int getMaxCount();

    String getValueStringAt(int pos);

    boolean isMismatch(int pos, byte ref, String chr, float snpThreshold);

    BisulfiteCounts getBisulfiteCounts();

    void finish();

    /**
     * Return the result of merging this alignment with {@code other}.
     * This alignment is not changed
     *
     * @param other
     * @param bisulfiteContext
     * @return
     */
    AlignmentCounts merge(AlignmentCounts other, AlignmentTrack.BisulfiteContext bisulfiteContext);

    static interface PositionIterator {
        int nextPosition();
    }
}
