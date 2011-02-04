/**
 * Copyright (c) 2011 by Fred Hutchinson Cancer Research Center.  All Rights Reserved.

 * This software is licensed under the terms of the GNU Lesser General
 * Public License (LGPL), Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.

 * THE SOFTWARE IS PROVIDED "AS IS." FRED HUTCHINSON CANCER RESEARCH CENTER MAKES NO
 * REPRESENTATIONS OR WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED,
 * INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 * PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS,
 * WHETHER OR NOT DISCOVERABLE.  IN NO EVENT SHALL FRED HUTCHINSON CANCER RESEARCH
 * CENTER OR ITS TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR
 * ANY DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR
 * CONSEQUENTIAL DAMAGES, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS,
 * REGARDLESS OF  WHETHER FRED HUTCHINSON CANCER RESEARCH CENTER SHALL BE ADVISED,
 * SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE
 * FOREGOING.
 */

package org.broad.igv.track;

import org.apache.log4j.Logger;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.renderer.SpliceJunctionRenderer;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.tribble.Feature;

import java.util.*;

/**
 * @author dhmay
 * @date Feb 3, 2011
 *
 * This class is a subclass of PackedFeatures that is for display of splice junctions. It overrides some
 * methods in order to deviate from superclass in two ways:
 * 1.  Features are allowed to be on the same line if flanking regions overlap
 * 2.  Features are ordered from top to bottom in ascending order of read depth
 * The goal is for the dominant isoform to appear on the top row.
 *
 * I think I've got the 90% case covered, but there's quite a bit of ambiguity here.  Some items for future work:
 * -The feature ordering will probably fail in certain conditions, e.g., an exon removal in which the situation
 * where the exon is present has more coverage in the first junction, but the absent-exon condition has more
 * coverage overall.  These are kind of degenerate cases, so only worth handling if someone complains.
 * -+ and - strand features occurring simultaneously.  Currently, not attempt is made to separate those features.
 * This could be done by dividing them up first, doing the packFeatures thing, and recombining.  This should probably
 * be done, but it's work.
 */
public class PackedFeaturesSpliceJunctions<T extends Feature> extends PackedFeatures{
    private static Logger log = Logger.getLogger(PackedFeaturesSpliceJunctions.class);

    PackedFeaturesSpliceJunctions(String chr, int start, int end, Iterator<T> iter, String trackName) {
        super(chr, start, end, iter, trackName);
    }

    /**
     * Splice junction features should be rendered on the same line even if their flanking regions overlap
     * @param feature
     * @return
     */
    protected int getFeatureStartForPacking(Feature feature)
    {
        return SpliceJunctionRenderer.getJunctionStart((IGVFeature) feature);
    }


    /**
     * Splice junction features should be rendered on the same line even if their flanking regions overlap
     * @param feature
     * @return
     */
    protected int getFeatureEndForPacking(Feature feature)
    {
        return SpliceJunctionRenderer.getJunctionEnd((IGVFeature) feature);
    }

    int getRowCount() {
        return getRows().size();
    }

    /**
     * Allocates each alignment to the rows such that there is no overlap. For splice junctions, priority queues
     * are ordered by feature score (read depth).  For the superclass, this is done by length
     *
     * @param iter TabixLineReader wrapping the collection of alignments
     */
    List<FeatureRow> packFeatures(Iterator iter) {

        List<FeatureRow> rows = new ArrayList(10);
        if (iter == null || !iter.hasNext()) {
            return rows;
        }

        maxFeatureLength = 0;
        int totalCount = 0;

        LinkedHashMap<Integer, PriorityQueue<T>> bucketArray = new LinkedHashMap();
        Comparator pqComparator = new Comparator<BasicFeature>() {
            public int compare(BasicFeature row1, BasicFeature row2) {
                return (int) (((IGVFeature) row2).getScore() - ((IGVFeature) row1).getScore());
            }
        };

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

        // Allocate features to rows
        FeatureRow currentRow = new FeatureRow();
        int allocatedCount = 1;
        int nextStart = currentRow.end + FeatureTrack.MINIMUM_FEATURE_SPACING;


        int lastKey = 0;
        int lastAllocatedCount = 0;
        while (allocatedCount < totalCount && rows.size() < maxLevels) {

            // Check to prevent infinite loops
            if (lastAllocatedCount == allocatedCount) {
                String msg = "Infinite loop detected while packing features for track: " + getTrackName() +
                        ".<br>Not all features will be shown." +
                        "<br>Please contact igv-help@broadinstitute.org";

                log.error(msg);
                MessageUtils.showMessage(msg);
                break;
            }
            lastAllocatedCount = allocatedCount;

            // Loop through alignments until we reach the end of the interval

            PriorityQueue<T> bucket = null;

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
}
