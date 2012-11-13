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

package org.broad.igv.track;

import org.apache.log4j.Logger;
import org.broad.igv.data.Interval;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.tribble.Feature;

import java.util.*;

/**
 * Represents a table of features, packed so there is no overlap.
 * Features are packed into rows, accessible via {@link #getRows}
 *
 * @author jrobinso
 * @date Oct 7, 2010
 */
public class PackedFeatures<T extends Feature>{
    protected String trackName;
    protected String chr;
    protected int start;
    protected int end;
    protected List<T> features;
    protected List<FeatureRow> rows;
    private static Logger log = Logger.getLogger(PackedFeatures.class);
    protected int maxFeatureLength = 0;
    protected static int maxLevels = 200;

    /**
     * No-arg constructor to allow subclassing
     */
    PackedFeatures(){
    }

    PackedFeatures(String chr, int start, int end) {
        this.chr = chr;
        this.start = start;
        this.end = end;
        features = Collections.emptyList();
        rows = Collections.emptyList();
    }

    PackedFeatures(String chr, int start, int end, Iterator<T> iter, String trackName) {
        this.trackName = trackName;
        this.chr = chr;
        this.start = start;
        this.end = end;
        features = new ArrayList(1000);
        rows = packFeatures(iter);
    }


    /**
     * Some types of Features (splice junctions) should be packed on the same row even if start and end overlap.
     * This can be overridden in a subclass
     * @param feature
     * @return
     */
    protected int getFeatureStartForPacking(Feature feature)
    {
        return feature.getStart();
    }


    /**
     * Some types of Features (splice junctions) should be packed on the same row even if start and end overlap.
     * This can be overridden in a subclass
     * @param feature
     * @return
     */
    protected int getFeatureEndForPacking(Feature feature)
    {
        return feature.getEnd();
    }

    int getRowCount() {
        return getRows().size();
    }

    public boolean containsInterval(String chr, int start, int end) {
        return this.getChr().equals(chr) && start >= this.getStart() && end <= this.getEnd();
    }

    /**
     * Allocates each feature to the rows such that there is no overlap.
     *
     * @param iter TabixLineReader wrapping the collection of alignments. Note that this should
     * really be an Iterator<T>, but it can't be subclassed if that's the case.
     */
    List<FeatureRow> packFeatures(Iterator iter) {

        List<FeatureRow> rows = new ArrayList(10);
        if (iter == null || !iter.hasNext()) {
            return rows;
        }

        maxFeatureLength = 0;
        int totalCount = 0;

        LinkedHashMap<Integer, PriorityQueue<T>> bucketArray = new LinkedHashMap();
        Comparator pqComparator = new Comparator<T>() {
            public int compare(Feature row1, Feature row2) {
                return (row2.getEnd() - row2.getStart()) - (row1.getEnd() - row2.getStart());
            }
        };

        // Allocate features to buckets,  1 bucket per base position
        while (iter.hasNext()) {
            T feature = (T) iter.next();
            maxFeatureLength = Math.max(maxFeatureLength,
                    getFeatureEndForPacking(feature) - getFeatureStartForPacking(feature));
            features.add(feature);

            int bucketNumber = getFeatureStartForPacking(feature);

            PriorityQueue<T> bucket = bucketArray.get(bucketNumber);
            if (bucket == null) {
                bucket = new PriorityQueue<T>(5, pqComparator);
                bucketArray.put(bucketNumber, bucket);
            }
            bucket.add(feature);
            totalCount++;

        }

        // Allocate features to rows, pulling at most 1 per bucket for each row
        FeatureRow currentRow = new FeatureRow();
        int allocatedCount = 0;
        int nextStart = Integer.MIN_VALUE;


        int lastKey = 0;
        int lastAllocatedCount = -1;
        while (allocatedCount < totalCount && rows.size() < maxLevels) {

            // Check to prevent infinite loops
            if (lastAllocatedCount == allocatedCount) {

                if(IGV.hasInstance()) {
                    String msg = "Infinite loop detected while packing features for track: " + getTrackName() +
                            ".<br>Not all features will be shown." +
                            "<br>Please contact igv-team@broadinstitute.org";

                    log.error(msg);
                    MessageUtils.showMessage(msg);
                }
                break;
            }
            lastAllocatedCount = allocatedCount;

            // Next row Loop through alignments until we reach the end of the interval

            PriorityQueue<T> bucket = null;
            // Advance to nextLine occupied bucket

            ArrayList<Integer> emptyBucketKeys = new ArrayList();
            for (Integer key : bucketArray.keySet()) {
                //if (key < lastKey) {
                //    String msg = "Features from track: " + trackName + " are not sorted.  Some features might not be shown.<br>" +
                //            "Please notify igv-help@broadinstitute.org";
                //    MessageUtils.showMessage(msg);
                //}
                lastKey = key;
                if (key >= nextStart) {
                    bucket = bucketArray.get(key);

                    T feature = bucket.poll();

                    if (bucket.isEmpty()) {
                        emptyBucketKeys.add(key);
                    }
                    currentRow.addFeature(feature);
                    nextStart = currentRow.end + FeatureTrack.MINIMUM_FEATURE_SPACING;
                    allocatedCount++;
                }
            }
            for (Integer key : emptyBucketKeys) {
                bucketArray.remove(key);
            }


            // We've reached the end of the interval,  start a new row
            if (currentRow.features.size() > 0) {
                rows.add(currentRow);
                lastAllocatedCount = 0;
            }
            currentRow = new FeatureRow();
            nextStart = 0;
            lastKey = 0;


        }
        // Add the last row
        if (currentRow.features.size() > 0) {
            rows.add(currentRow);
        }

        return rows;
    }

    public String getTrackName() {
        return trackName;
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

    public List<T> getFeatures() {
        return features;
    }

    public List<FeatureRow> getRows() {
        return rows;
    }

    public int getMaxFeatureLength() {
        return maxFeatureLength;
    }

    class FeatureRow {
        int start;
        int end;
        List<T> features;

        public FeatureRow() {
            this.features = new ArrayList(100);
        }

        public void addFeature(T feature) {
            if (features.isEmpty()) {
                this.start = getFeatureStartForPacking(feature);
            }
            features.add(feature);
            end = getFeatureEndForPacking(feature);
        }

        public List<T> getFeatures() {
            return features;
        }
    }
}
