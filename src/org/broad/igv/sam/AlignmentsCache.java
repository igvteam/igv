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

import org.broad.igv.feature.Range;

import java.util.*;

/**
 * Stores {@link AlignmentInterval}s for easy lookup by position
 *
 * @author jacob, jrobinso
 * @date 2014-Jan-10
 */
public class AlignmentsCache {

    Map<String, List<AlignmentInterval>> intervals;

    AlignmentsCache() {
        intervals = new HashMap<String, List<AlignmentInterval>>();
    }

    AlignmentsCache(AlignmentsCache cache){
        this.intervals = new HashMap<String, List<AlignmentInterval>>(cache.intervals);
    }

    /** Add the specified interval to the cache. Replaces any existing interval
     * which it fully contains.
     * @param interval
     * @return The old interval, null if it didn't exist
     */
    AlignmentInterval put(AlignmentInterval interval) {
        String chr = interval.getChr();
        List<AlignmentInterval> iList = intervals.get(chr);
        int currentIndex = -1;
        if (iList == null) {
            iList = new ArrayList<AlignmentInterval>();
            intervals.put(chr, iList);
        } else {
            currentIndex = getIndexOf(chr, interval.getStart(), interval.getEnd());
        }

        if (currentIndex >= 0) {
            AlignmentInterval currentInterval = iList.get(currentIndex);
            iList.set(currentIndex, interval);
            return currentInterval;
        } else {
            iList.add(interval);
            return null;
        }

    }

    AlignmentInterval get(Range range) {
        String chr = range.getChr();
        int index = getIndexOf(chr, range.getStart(), range.getEnd());
        return index >= 0 ? intervals.get(chr).get(index) : null;
    }

    AlignmentInterval get(String chr, int start, int end) {
        int index = getIndexOf(chr, start, end);
        return index >= 0 ? intervals.get(chr).get(index) : null;
    }

    private int getIndexOf(String chr, int start, int end) {

        List<AlignmentInterval> iList = intervals.get(chr);
        if (iList != null) {
            for (int ii = 0; ii < iList.size(); ii++) {
                AlignmentInterval ai = iList.get(ii);
                if (ai.contains(chr, start, end)) {
                    return ii;
                }
            }
        }
        return -1;

    }

    public Collection<AlignmentInterval> getLoadedIntervals() {
        List<AlignmentInterval> alignmentIntervals = new ArrayList<AlignmentInterval>();
        for(List<AlignmentInterval> intervals: this.intervals.values()){
            alignmentIntervals.addAll(intervals);
        }
        return alignmentIntervals;
    }

    public void clear() {
        this.intervals.clear();
    }
}
