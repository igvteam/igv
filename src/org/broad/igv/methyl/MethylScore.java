/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not
 * responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL), Version 2.1 which is
 * available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.methyl;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.Strand;
import org.broad.igv.track.WindowFunction;

/**
 * @author Jim Robinson
 * @date 4/19/12
 */
public class MethylScore implements LocusScore {

    String chr;
    int start;
    int end;
    Strand strand;
    float percentMethylated;
    int totalCount;

    public MethylScore(String chr, int start, int end, Strand strand, float percentMethylated, int totalCount) {
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.strand = strand;
        this.percentMethylated = percentMethylated;
        this.totalCount = totalCount;
    }

    public String getChr() {
        return chr;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public Strand getStrand() {
        return strand;
    }

    public float getScore() {
        return percentMethylated;
    }

    public int getCount() {
        return totalCount;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    /**
     * Return a string to be used for popup text.   The WindowFunction is passed
     * in so it can be used t annotate the value.  The LocusScore object itself
     * does not "know" from what window function it was derived
     *
     * @param position       Zero-based genome position
     * @param windowFunction
     * @return
     */
    public String getValueString(double position, WindowFunction windowFunction) {
        return percentMethylated + "%" + " [" + totalCount + "]" +
                (strand == Strand.POSITIVE ? " (+)" : (strand == Strand.NEGATIVE ? " (-)" : ""));
    }

}
