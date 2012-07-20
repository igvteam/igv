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

package org.broad.igv.util.collections;

import org.apache.commons.collections.Predicate;
import org.broad.igv.data.Interval;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.util.Utilities;

import java.util.*;

/**
 * Used for caching data oriented in a genomic interval.
 * Data is stored by sequence name (string), start, end, and zoom.
 * <p/>
 * Lookups are not designed to be efficient (not a B-tree)
 * User: jacob
 * Date: 2012-Jul-13
 */
public class CachedIntervals<T extends Interval> {

    private Map<String, List<T>> map = Collections.synchronizedMap(new HashMap<String, List<T>>());

    private static final int DEFAULT_CACHE_SIZE = 5;
    private static final int DEFAULT_MAX_INTERVAL_SIZE = (int) 10e6;

    protected int cacheSize = DEFAULT_CACHE_SIZE;
    protected int maxIntervalSize = DEFAULT_MAX_INTERVAL_SIZE;

    public CachedIntervals() {
    }

    public CachedIntervals(int cacheSize, int maxIntervalSize) {
        this.cacheSize = cacheSize;
        this.maxIntervalSize = maxIntervalSize;
    }

    /**
     * Get the list of intervals which contain the specified interval,
     * sorted by start position.
     *
     * @param chr
     * @param start
     * @param end
     * @param zoom
     * @return
     */
    public List<T> getContains(final String chr, final int start, final int end, final int zoom) {
        Predicate<Interval> pred = new Predicate<Interval>() {
            @Override
            public boolean evaluate(Interval interval) {
                return interval.contains(chr, start, end, zoom);
            }
        };
        return getGen(chr, pred);
    }

    /**
     * Get the list of Intervals which overlap the specified interval, sorted by start position.
     *
     * @param chr
     * @param start
     * @param end
     * @param zoom
     * @return
     */
    public List<T> getOverlaps(final String chr, final int start, final int end, final int zoom) {
        Predicate<Interval> pred = new Predicate<Interval>() {
            @Override
            public boolean evaluate(Interval interval) {
                return interval.overlaps(chr, start, end, zoom);
            }
        };
        return getGen(chr, pred);
    }

    /**
     * Return all stored intervals for the specified seqName, sorted
     * by start position.
     *
     * @param seqName
     * @return
     */
    public List<T> get(String seqName) {
        return getGen(seqName, null);
    }


    /**
     * Retrieve stored intervals based which match predicate, sorted by start
     * If predicate is null, return all intervals for the specified chromosome
     *
     * @param seqName
     * @param predicate
     * @return
     */
    private List<T> getGen(String seqName, Predicate<Interval> predicate) {
        List<T> intervals = map.get(seqName);

        if (intervals == null) return null;

        List<T> returnedIntervals = new ArrayList<T>(intervals);
        if (predicate != null) {
            Utilities.filter(returnedIntervals, predicate);
        }
        FeatureUtils.sortFeatureList(returnedIntervals);

        return returnedIntervals;
    }

    /**
     * Stores this interval. addingInterval.getChr() is used
     * as the key. This interval will be available in
     * get(addingInterval.getChr()), although it may have
     * been merged into an existing interval.
     * <p/>
     * This method is synchronized because the internal contents
     * of the map may be modified as a result, if intervals
     * are merged.
     *
     * @param addingInterval
     */
    public synchronized void put(T addingInterval) {
        boolean doAdd = true;
        String key = addingInterval.getChr();
        List<T> intervals = map.get(key);
        if (intervals == null) {
            intervals = new LinkedList<T>();
            map.put(key, intervals);
        } else {
            List<T> overlaps = getOverlaps(addingInterval.getChr(), addingInterval.getStart(), addingInterval.getEnd(),
                    addingInterval.getZoom());

            //If it overlaps an existing interval, merge it in.
            //Otherwise, just add it to the end
            if (overlaps != null && overlaps.size() > 0) {
                Interval overlap = overlaps.get(0);
                overlap.merge(addingInterval);
                doAdd = false;

                //Prevent interval from growing without bound
                int intervalSize = overlap.getEnd() - overlap.getStart();
                if (intervalSize > maxIntervalSize) {
                    int newStart = Math.min(overlap.getStart(), addingInterval.getStart());
                    int newEnd = Math.max(overlap.getEnd(), overlap.getEnd());
                    overlap.trimTo(overlap.getChr(), newStart, newEnd, overlap.getZoom());
                }
            }
        }

        if (doAdd)
            intervals.add(addingInterval);


        while (intervals.size() > cacheSize) {
            intervals.remove(0);
        }


    }


    public int getMaxIntervalSize() {
        return maxIntervalSize;
    }

    public void setMaxIntervalSize(int maxIntervalSize) {
        this.maxIntervalSize = maxIntervalSize < 0 ? Integer.MAX_VALUE : maxIntervalSize;
    }

    public int size() {
        return map.size();
    }

    public Collection<T> getLoadedIntervals() {
        List<T> allLoadedIntervals = new ArrayList<T>(map.size());
        synchronized (map) {
            for (Map.Entry<String, List<T>> entry : map.entrySet()) {
                allLoadedIntervals.addAll(entry.getValue());
            }
        }
        return allLoadedIntervals;
    }

    public void clear() {
        map.clear();
    }
}
