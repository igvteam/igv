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

package org.broad.igv.sam;

import org.broad.igv.feature.Range;
import org.broad.igv.util.collections.LRUCache;

import java.util.ArrayList;
import java.util.Collection;

/**
 * Stores objects by position so they can be looked up by inexact position (contains)
 *
 * @author jacob, jrobinso
 * @date 2014-Jan-10
 */
class PositionCache<V> {

    private LRUCache<Range, V> intervals;

    private static final int MIN_MAX_ENTRIES = 10;

    public PositionCache() {
        intervals = new LRUCache<Range, V>(MIN_MAX_ENTRIES);
    }

    public PositionCache(PositionCache<V> cache){
        this.intervals = new LRUCache<Range, V>(MIN_MAX_ENTRIES);
        this.intervals.putAll(cache.intervals);
    }

    /** Add the specified interval to the cache. Replaces any existing interval
     * which contains the given range
     * @param range
     * @param value
     * @return The old interval, null if it didn't exist
     */
    public V put(Range range, V value) {
        Range currentRangeKey = getKeyForRange(range);
        Range keyToUse = currentRangeKey != null ? currentRangeKey : range;
        return intervals.put(keyToUse, value);

    }

    public V getForRange(Range range){
        Range key = getKeyForRange(range);
        return key != null ? intervals.get(key) : null;
    }

    private Range getKeyForRange(Range range) {
        String chr = range.getChr();
        for (Range cachedRange : intervals.keySet()) {
            if (cachedRange.contains(chr, range.getStart(), range.getEnd())) {
                return cachedRange;
            }
        }
        return null;
    }


    public boolean containsRange(Range range) {
        return getKeyForRange(range) != null;
    }

    public Collection<V> values() {
        return new ArrayList<V>(this.intervals.values());
    }

    public void clear() {
        this.intervals.clear();
    }

    public synchronized void setMaxEntries(int inMaxEntries){
        int newMax = Math.max(MIN_MAX_ENTRIES, inMaxEntries);
        intervals.setMaxEntries(newMax);
    }

}
