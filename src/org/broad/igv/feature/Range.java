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

package org.broad.igv.feature;

import com.google.java.contract.util.Objects;
import org.broad.tribble.Feature;


/**
 * Basic class to specify a genomic interval.
 * Coordinates are intended to be 0-based half-open
 * Please do not add or remove any fields (want to keep it very simple).
 * Additional methods for calculating overlap are okay
 *
 * @author jacob
 * @author 2013-May-20
 */
public class Range implements Feature {

    protected String chr = null;
    protected int start = -1;
    protected int end = -1;

    public Range(String chr, int start, int end){
        this.chr = chr;
        this.start = start;
        this.end = end;
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

    public int getLength(){
        return end - start;
    }

    /**
     * Determine whether this interval fully contains the specified
     * input interval.
     *
     * A negative input start position has special meaning.  It is considered within the interval if the interval
     * contains position "0".
     *
     * @param chr
     * @param start
     * @param end
     * @return
     */
    public boolean contains(String chr, int start, int end) {
        return Objects.equal(this.chr, chr)
                && this.start <= (start < 0 ? 0 : start)
                && this.end >= end;
    }

    /**
     * Determine whether there is any overlap between this interval and the specified interval
     *
     * Negative positions have no special meaning
     * @param chr
     * @param start
     * @param end
     * @return
     */
    public boolean overlaps(String chr, int start, int end) {
        return Objects.equal(this.chr, chr) && this.start <= end && this.end >= start;
    }

    public boolean overlaps(Range range) {
        return this.overlaps(range.getChr(), range.getStart(), range.getEnd());
    }

}
