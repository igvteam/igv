/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
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
 * Genomic interval showing which areas have had reads removed (downsampled)
* @author jrobinso
*         Date: 8/16/12
*         Time: 10:03 PM
*/
public class DownsampledInterval implements Feature {
    private int start;
    private int end;
    private int count;

    public DownsampledInterval(int start, int end, int count) {
        this.start = start;
        this.end = end;
        this.count = count;
    }

    public String toString() {
        return start + "-" + end + " (" + count + ")";
    }

    public int getCount() {
        return count;
    }

    public int getEnd() {
        return end;
    }

    public int getStart() {
        return start;
    }

    public String getChr() {
        return null;
    }

    public String getValueString() {
        return "Interval [" + start + "-" + end + "] <br>" + count + " reads removed.";
    }

    public void incCount() {
        this.count += 1;
    }
}
