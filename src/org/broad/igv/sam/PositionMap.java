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
 * Stores objects by position so they can be looked up by inexact position (contains)
 *
 * @author jacob, jrobinso
 * @date 2014-Jan-10
 */
class PositionMap<V> {

    /**
     * Map from chr -> objects
     */
    Map<String, List<ValuedRange>> intervals;

    public PositionMap() {
        intervals = new HashMap<String, List<ValuedRange>>();
    }

    public PositionMap(PositionMap<V> cache){
        this.intervals = new HashMap<String, List<ValuedRange>>(cache.intervals);
    }

    /** Add the specified interval to the cache. Replaces any existing interval
     * which it fully contains.
     * @param range
     * @param value
     * @return The old interval, null if it didn't exist
     */
    public V put(Range range, V value) {
        String chr = range.getChr();
        List<ValuedRange> iList = intervals.get(chr);
        int currentIndex = -1;
        if (iList == null) {
            iList = new ArrayList<ValuedRange>();
            intervals.put(chr, iList);
        } else {
            currentIndex = getIndexOf(range);
        }

        ValuedRange newValue = new ValuedRange(range, value);
        if (currentIndex >= 0) {
            ValuedRange currentValue = iList.get(currentIndex);
            iList.set(currentIndex, newValue);
            return currentValue.value;
        } else {
            iList.add(newValue);
            return null;
        }

    }

    public V get(Range range) {
        String chr = range.getChr();
        int index = getIndexOf(range);
        return index >= 0 ? intervals.get(chr).get(index).value : null;
    }

    private int getIndexOf(Range range) {
        String chr = range.getChr();
        List<ValuedRange> iList = intervals.get(chr);
        if (iList != null) {
            for (int ii = 0; ii < iList.size(); ii++) {
                ValuedRange vr = iList.get(ii);
                if (vr.contains(chr, range.getStart(), range.getEnd())) {
                    return ii;
                }
            }
        }
        return -1;
    }


    public boolean contains(Range range) {
        return getIndexOf(range) >= 0;
    }

    public Collection<V> values() {
        List<V> values = new ArrayList<V>();
        for(List<ValuedRange> interval: this.intervals.values()){
            for(ValuedRange vr: interval){
                values.add(vr.value);
            }
        }
        return values;
    }

    public void clear() {
        this.intervals.clear();
    }

    private class ValuedRange extends Range{

        private V value;

        ValuedRange(Range range, V value){
            super(range.getChr(), range.getStart(), range.getEnd());
            this.value = value;
        }
    }
}
