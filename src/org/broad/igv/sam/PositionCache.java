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
import org.broad.igv.util.collections.LRUCache;

import java.util.*;

/**
 * Stores objects by position so they can be looked up by inexact position (contains)
 *
 * @author jacob, jrobinso
 * @date 2014-Jan-10
 */
class PositionCache<V> {

    /**
     * Map from chr -> caches
     */
    private Map<String, LRUCache<Range, V>> intervals;

    private static final int MIN_MAX_ENTRIES = 4;
    private int maxEntriesPerChr = MIN_MAX_ENTRIES;

    public PositionCache() {
        intervals = Collections.synchronizedMap(new HashMap<String, LRUCache<Range, V>>());
    }

    public PositionCache(PositionCache<V> cache){
        this.intervals = Collections.synchronizedMap(new HashMap<String, LRUCache<Range, V>>(cache.intervals));
    }

    /** Add the specified interval to the cache. Replaces any existing interval
     * which it fully contains.
     * @param range
     * @param value
     * @return The old interval, null if it didn't exist
     */
    public synchronized V put(Range range, V value) {
        String chr = range.getChr();
        LRUCache<Range, V> chrVs = intervals.get(chr);
        Range currentRangeKey = null;
        if (chrVs == null) {
            chrVs = new LRUCache<Range, V>(this.maxEntriesPerChr);
            intervals.put(chr, chrVs);
        } else {
            currentRangeKey = getKeyFor(range);
        }

        if (currentRangeKey != null) {
            return chrVs.put(currentRangeKey, value);
        } else {
            chrVs.put(range, value);
            return null;
        }

    }

    public V get(Range range) {
        String chr = range.getChr();
        Range keyRange = getKeyFor(range);
        return keyRange != null ? intervals.get(chr).get(keyRange) : null;
    }

    private Range getKeyFor(Range range) {
        String chr = range.getChr();
        LRUCache<Range, V> chrVs = intervals.get(chr);
        if (chrVs != null) {
            for (Range cachedRange: chrVs.keySet()) {
                if (cachedRange.contains(chr, range.getStart(), range.getEnd())) {
                    return cachedRange;
                }
            }
        }
        return null;
    }


    public boolean contains(Range range) {
        return getKeyFor(range) != null;
    }

    public Collection<V> values() {
        List<V> values = new ArrayList<V>();
        for(LRUCache<Range, V> interval: this.intervals.values()){
            values.addAll(interval.values());
        }
        return values;
    }

    public void clear() {
        this.intervals.clear();
    }

    public synchronized void setMaxEntriesPerChr(int inMaxEntries){
        this.maxEntriesPerChr = Math.max(MIN_MAX_ENTRIES, inMaxEntries);
        for(LRUCache<Range, V> cache: this.intervals.values()){
            cache.setMaxEntries(this.maxEntriesPerChr);
        }
    }

}
