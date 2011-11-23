/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam;

import org.apache.log4j.Logger;
import org.broad.igv.sam.AlignmentInterval.Row;

import javax.sound.midi.SysexMessage;
import java.util.*;

/**
 * A utility class to experiment with alignment packing methods.
 * <p/>
 * Numbers:
 * <p/>
 * packAlignments1  (original)
 * Packed  19075 out of 19075 in 59 rows in:     0.046 seconds
 * Packed  111748 out of 431053 in 1000 rows in: 14.313 seconds
 * <p/>
 * packAlignments2  (row by row)
 * Packed 17027 out of 17027 in 58 rows:        0.035 seconds
 * Packed 113061 out of 431053 in 1000 rows in  8.276 seconds
 * :
 * packAlignments2b  (row by row with hash)
 * Packed 15274 out of 15274 in 58 rows in:     0.011 seconds
 * Packed 101595 out of 423736 in 1000 rows in: 0.177 seconds
 * <p/>
 * packAlignments3  (priority queue)
 * Packed 19075 out of 19075 in 63 rows in:      0.044 seconds
 * Packed 104251 out of 430716 in 1000 rows in:  0.108 seconds
 *
 * @author jrobinso
 */
public class AlignmentPacker {

    private static Logger log = Logger.getLogger(AlignmentPacker.class);

    /**
     * Minimum gap between the end of one alignment and start of another.
     */
    public static final int MIN_ALIGNMENT_SPACING = 5;
    private static final int MAX_ROWS = 100000;
    private Comparator lengthComparator;

    public AlignmentPacker() {
        lengthComparator = new Comparator<Alignment>() {
            public int compare(Alignment row1, Alignment row2) {
                return (row2.getEnd() - row2.getStart()) -
                        (row1.getEnd() - row2.getStart());

            }
        };
    }


    /**
     * Allocates each alignment to the rows such that there is no overlap.
     *
     * @param iter Iterator wrapping the collection of alignments
     */
    public List<AlignmentInterval.Row> packAlignments(
            Iterator<Alignment> iter,
            int end,
            boolean pairAlignments,
            AlignmentTrack.SortOption groupBy,
            int maxLevels) {


        List<Row> alignmentRows = new ArrayList(1000);
        if (iter == null || !iter.hasNext()) {
            return alignmentRows;
        }

        if (groupBy == null) {
            pack(iter, end, pairAlignments, lengthComparator, alignmentRows, maxLevels);
        } else {
            // Separate by group
            List<Alignment> nullGroup = new ArrayList();
            HashMap<String, List<Alignment>> groupedAlignments = new HashMap();
            while (iter.hasNext()) {
                Alignment al = iter.next();
                String groupKey = getGroupValue(al, groupBy);
                if (groupKey == null) nullGroup.add(al);
                else {
                    List<Alignment> group = groupedAlignments.get(groupKey);
                    if (group == null) {
                        group = new ArrayList(1000);
                        groupedAlignments.put(groupKey, group);
                    }
                    group.add(al);
                }
            }
            List<String> keys = new ArrayList(groupedAlignments.keySet());
            Collections.sort(keys);
            for (String key : keys) {
                List<Alignment> group = groupedAlignments.get(key);
                pack(group.iterator(), end, pairAlignments, lengthComparator, alignmentRows, maxLevels);
            }
            pack(nullGroup.iterator(), end, pairAlignments, lengthComparator, alignmentRows, maxLevels);
        }

        return alignmentRows;

    }

    private String getGroupValue(Alignment al, AlignmentTrack.SortOption groupBy) {
        switch (groupBy) {

            case STRAND:
                return String.valueOf(al.isNegativeStrand());
            case SAMPLE:
                return al.getSample();
            case READ_GROUP:
                return al.getReadGroup();
        }
        return null;
    }

    private void pack(Iterator<Alignment> iter, int end, boolean pairAlignments, Comparator lengthComparator, List<Row> alignmentRows,
                      int maxLevels) {

        if (!iter.hasNext()) {
            return;
        }

        Map<String, PairedAlignment> pairs = null;
        if (pairAlignments) {
            pairs = new HashMap(1000);
        }


        // Strictly speaking we should loop discarding dupes, etc.
        Alignment firstAlignment = iter.next();
        if (pairAlignments && firstAlignment.isPaired() && firstAlignment.isProperPair() && firstAlignment.getMate().isMapped()) {
            String readName = firstAlignment.getReadName();
            PairedAlignment pair = new PairedAlignment(firstAlignment);
            pairs.put(readName, pair);
            firstAlignment = pair;
        }

        int start = firstAlignment.getStart();
        int bucketCount = end - start + 1;

        // Create buckets.  We use priority queues to keep the buckets sorted by alignment length.  However this
        // is probably a neeedless complication,  any collection type would do.
        PriorityQueue firstBucket = new PriorityQueue(5, lengthComparator);
        firstBucket.add(firstAlignment);

        // Use dense buckets for < 50,000, sparse otherwise

        BucketCollection buckets;
        if (bucketCount < 50000) {
            buckets = new DenseBucketCollection(bucketCount);
        } else {
            buckets = new SparseBucketCollection();
        }
        buckets.set(0, firstBucket);

        int totalCount = 1;

        //  Allocate alignments to buckets based on position

        while (iter.hasNext()) {

            Alignment al = iter.next();
            String readName = al.getReadName();

            if (al.isMapped()) {

                Alignment alignment = al;
                if (pairAlignments && al.isPaired() && al.getMate().isMapped() && al.getChr().equals(al.getMate().getChr())) {

                    PairedAlignment pair = pairs.get(readName);
                    if (pair == null) {
                        pair = new PairedAlignment(al);
                        pairs.put(readName, pair);
                        alignment = pair;
                    } else {
                        if (al.getChr().equals(pair.getChr())) {
                            // Add second alignment to pair
                            pair.setSecondAlignment(al);
                            pairs.remove(readName);
                            continue;
                        }

                    }
                }


                // We can get negative buckets if softclipping is on as the alignments are only approximately
                // sorted.  Throw all alignments < start in the first bucket.
                int bucketNumber = Math.max(0, alignment.getStart() - start);
                if (bucketNumber < bucketCount) {
                    PriorityQueue bucket = buckets.get(bucketNumber);
                    if (bucket == null) {
                        bucket = new PriorityQueue<Alignment>(5, lengthComparator);
                        buckets.set(bucketNumber, bucket);
                    }
                    bucket.add(alignment);
                    totalCount++;
                } else {
                    log.debug("Alignment out of bounds: " + alignment.getStart() + " (> " + end);
                }


            }
        }

        buckets.finishedAdding();

        // Allocate alignments to rows
        long t0 = System.currentTimeMillis();
        int allocatedCount = 0;
        int nextStart = start;
        Row currentRow = new Row();
        while (allocatedCount < totalCount) { // && alignmentRows.size() < maxLevels) {

            // Loop through alignments until we reach the end of the interval
            while (nextStart <= end) {
                PriorityQueue<Alignment> bucket;

                // Advance to next occupied bucket
                int bucketNumber = nextStart - start;
                bucket = buckets.getNextBucket(bucketNumber);

                // Pull the next alignment out of the bucket and add to the current row
                if (bucket == null) {
                    break;
                } else {
                    Alignment alignment = bucket.remove();
                    currentRow.addAlignment(alignment);
                    nextStart = currentRow.getLastEnd() + MIN_ALIGNMENT_SPACING;
                    allocatedCount++;
                }

            }

            // We've reached the end of the interval,  start a new row
            if (currentRow.alignments.size() > 0) {
                alignmentRows.add(currentRow);
            }

            //if (alignmentRows.size() >= maxLevels) {
            //    currentRow = null;
            //    break;
            //}

            currentRow = new Row();
            nextStart = start;
        }
        if (log.isDebugEnabled()) {
            long dt = System.currentTimeMillis() - t0;
            log.debug("Packed alignments in " + dt);
        }

        // Add the last row
        if (currentRow != null && currentRow.alignments.size() > 0) {
            alignmentRows.add(currentRow);
        }

    }


    static interface BucketCollection {

        void set(int idx, PriorityQueue<Alignment> bucket);

        PriorityQueue<Alignment> get(int idx);

        PriorityQueue<Alignment> getNextBucket(int bucketNumber);

        void finishedAdding();

    }


    static class DenseBucketCollection implements BucketCollection {

        int lastBucketNumber = -1;

        PriorityQueue<Alignment>[] bucketArray;

        DenseBucketCollection(int bucketCount) {
            bucketArray = new PriorityQueue[bucketCount];
        }

        public void set(int idx, PriorityQueue<Alignment> bucket) {
            bucketArray[idx] = bucket;
        }

        public PriorityQueue<Alignment> get(int idx) {
            return bucketArray[idx];
        }


        /**
         * Return the next occupied bucket after bucketNumber
         *
         * @param bucketNumber
         * @return
         */
        public PriorityQueue<Alignment> getNextBucket(int bucketNumber) {

            if (bucketNumber == lastBucketNumber) {
                // TODO -- detect inf loop here
            }

            PriorityQueue<Alignment> bucket = null;
            while (bucketNumber < bucketArray.length) {
                bucket = bucketArray[bucketNumber];
                if (bucket != null) {
                    if (bucket.isEmpty()) {
                        bucketArray[bucketNumber] = null;
                    } else {
                        return bucket;
                    }
                }
                bucketNumber++;
            }
            return bucket;
        }

        public void finishedAdding() {
            // nothing to do
        }
    }


    static class SparseBucketCollection implements BucketCollection {

        boolean finished = false;
        List<Integer> keys;
        HashMap<Integer, PriorityQueue<Alignment>> buckets;

        SparseBucketCollection() {
            buckets = new HashMap(1000);
        }

        public void set(int idx, PriorityQueue<Alignment> bucket) {
            if (finished) {
                System.err.println("Error: bucket added after finishAdding() called");
            }
            buckets.put(idx, bucket);
        }

        public PriorityQueue<Alignment> get(int idx) {
            return buckets.get(idx);
        }

        /**
         * Return the next occupied bucket after bucketNumber
         *
         * @param bucketNumber
         * @return
         */
        public PriorityQueue<Alignment> getNextBucket(int bucketNumber) {


            int min = 0;
            int max = keys.size() - 1;
            int minBucket = keys.get(0);
            int maxBucket = keys.get(max);
            if (bucketNumber < minBucket) {
                buckets.get(minBucket);
            }
            if (bucketNumber > maxBucket) {
                return null;
            }

            int emptyCount = 0;
            while ((max - min) > 5) {
                int mid = (max + min) / 2;
                int key = keys.get(mid);

                // Skip empty buckets
                PriorityQueue<Alignment> bucket = buckets.get(key);
                if (bucket.isEmpty()) {
                    emptyCount++;
                    max--;
                } else {
                    if (key == bucketNumber) {
                        return bucket;
                    }
                    if (key > bucketNumber) {
                        max = mid;
                    } else {
                        min = mid;
                    }
                }
            }

            // Now march from min to max until we cross bucketNumber
            for (int i = min; i <buckets.size(); i++) {
                int key = keys.get(i);
                PriorityQueue<Alignment> bucket = buckets.get(key);
                if (bucket.isEmpty()) continue;

                // System.out.println("    key = " + key);
                if (key >= bucketNumber) {
                    //   System.out.println("        RETURN " + key);
                    //   System.out.println();
                    return buckets.get(key);
                }
            }

            // Cleanup if we have too many empty buckets
            if ((((double) emptyCount) / buckets.size()) > 0.3) {
                HashMap<Integer, PriorityQueue<Alignment>> newBuckets = new HashMap<Integer, PriorityQueue<Alignment>>(buckets.size());
                for (Map.Entry<Integer, PriorityQueue<Alignment>> entry : buckets.entrySet()) {
                    final PriorityQueue<Alignment> b = entry.getValue();
                    if (!b.isEmpty()) {
                        newBuckets.put(entry.getKey(), b);
                    }
                }
                buckets = newBuckets;
                finishedAdding();
            }

            return buckets.get(max);
        }

        public void finishedAdding() {
            finished = true;
            keys = new ArrayList<Integer>(buckets.keySet());
            Collections.sort(keys);
        }
    }

}
