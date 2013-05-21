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

package org.broad.igv.data;

import com.google.java.contract.util.Objects;
import org.broad.tribble.Feature;


/**
 * Basic class to specify a genomic interval.
 *
 * User: jacob
 * Date: 2013-May-20
 */
public class Interval implements Feature {

    protected String chr = null;
    protected int start = -1;
    protected int end = -1;

    public Interval(String chr, int start, int end){
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

    /**
     * Determine whether this interval fully contains the specified
     * input interval
     *
     * @param chr
     * @param start
     * @param end
     * @return
     */
    public boolean contains(String chr, int start, int end) {
        return Objects.equal(this.chr, chr) && this.start <= start && this.end >= end;
    }

    /**
     * Determine whether there is any overlap between this interval and the specified interval
     * @param chr
     * @param start
     * @param end
     * @return
     */
    public boolean overlaps(String chr, int start, int end) {
        return Objects.equal(this.chr, chr) && this.start <= end && this.end >= start;
    }



}
