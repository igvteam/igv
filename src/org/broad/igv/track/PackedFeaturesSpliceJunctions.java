/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Fred Hutchinson Cancer Research Center and Broad Institute
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


package org.broad.igv.track;

import org.apache.log4j.Logger;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.SpliceJunctionFeature;
import org.broad.igv.feature.Strand;
import org.broad.igv.renderer.SpliceJunctionRenderer;
import org.broad.igv.ui.util.MessageUtils;
import htsjdk.tribble.Feature;

import java.util.*;

/**
 * @author dhmay
 * @date Feb 3, 2011
 *
 * This class is a subclass of PackedFeatures that is for display of splice junctions. It overrides some
 * methods in order to deviate from superclass in two ways:
 * 1.  Features are allowed to be on the same line if flanking regions overlap
 * 2.  Overlapping features are allowed to be on the same line if they are from different strands
 * 3.  Features are ordered from top to bottom in ascending order of read depth
 * The goal is for the dominant isoform to appear on the top row.
 *
 * I think I've got the 90% case covered, but there's quite a bit of ambiguity here.  Some items for future work:
 * -The feature ordering will probably fail in certain conditions, e.g., an exon removal in which the situation
 * where the exon is present has more coverage in the first junction, but the absent-exon condition has more
 * coverage overall.  These are kind of degenerate cases, so only worth handling if someone complains.
 */
public class PackedFeaturesSpliceJunctions<T extends Feature> extends PackedFeatures {
    private static Logger log = Logger.getLogger(PackedFeaturesSpliceJunctions.class);

    public PackedFeaturesSpliceJunctions(String chr, int start, int end, Iterator<T> iter, String trackName) {
        super(chr, start, end, iter, trackName);
    }

    /**
     * Splice junction features should be rendered on the same line even if their flanking regions overlap
     * @param feature
     * @return
     */
    protected int getFeatureStartForPacking(Feature feature)
    {
        return ((SpliceJunctionFeature) feature).getJunctionStart();
    }


    /**
     * Splice junction features should be rendered on the same line even if their flanking regions overlap
     * @param feature
     * @return
     */
    protected int getFeatureEndForPacking(Feature feature)
    {
        return ((SpliceJunctionFeature) feature).getJunctionEnd();
    }

    int getRowCount() {
        return getRows().size();
    }

    /**
     * Allocates each alignment to the rows such that there is no overlap. For splice junctions, priority queues
     * are ordered by feature score (read depth).  For the superclass, this is done by length.
     * Since splice junctions only interfere with each other within a strand, break up the iterator into one
     * iterator per strand, farm out the work per strand, and reintegrate.
     *
     * This seems nice and clean, but it's actually not that efficient, since we're handling all of one strand
     * and then all of the other -- we're essentially buffering all the second strand's features until we get
     * to them
     *
     * @param iter TabixLineReader wrapping the collection of alignments
     */
    List<FeatureRow> packFeatures(Iterator iter) {
        IteratorSplitterByCharge iterSplitter = new IteratorSplitterByCharge(iter);

        List<FeatureRow> posRows = packFeaturesOneStrand(iterSplitter.getPosIter());
        List<FeatureRow> negativeRows = packFeaturesOneStrand(iterSplitter.getNegIter());

        Comparator startComparator = new Comparator<Feature>() {
            public int compare(Feature row1, Feature row2) {
                return row1.getStart() - row2.getStart();
            }
        };

        int numRows = Math.max(posRows.size(), negativeRows.size());
        List<FeatureRow> result = new ArrayList<FeatureRow>(numRows);
        features.clear();
        for (int i=0; i<numRows; i++)
        {
            List<Feature> posAndNegFeatures = new ArrayList<Feature>();
            if (negativeRows.size() > i)
                posAndNegFeatures.addAll(negativeRows.get(i).getFeatures());
            if (posRows.size() > i)
                posAndNegFeatures.addAll(posRows.get(i).getFeatures());

            if (!posAndNegFeatures.isEmpty())
            {
                Collections.sort(posAndNegFeatures, startComparator);
                FeatureRow resultRow = new FeatureRow();
                for (Feature feature : posAndNegFeatures)
                    resultRow.addFeature(feature);
                result.add(resultRow);
                features.addAll(posAndNegFeatures);
            }
        }
        Collections.sort(features, startComparator);
        return result;
    }


    /**
     * This does the real work of packing features, pretty much the same way the superclass does, except
     * that features can overlap in their flanking regions and they're ordered among the rows by score
     * @param iter
     * @return
     */
    List<FeatureRow> packFeaturesOneStrand(Iterator iter) {
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
        int allocatedCount = 0;
        int nextStart = currentRow.end + FeatureTrack.MINIMUM_FEATURE_SPACING;

        int lastAllocatedCount = -1;
        while (allocatedCount < totalCount && rows.size() < maxLevels) {

            // Check to prevent infinite loops
            if (lastAllocatedCount == allocatedCount) {
                String msg = "Infinite loop detected while packing features for track: " + getTrackName() +
                        ".<br>Not all features will be shown." +
                        "<br>Please contact igv-team@broadinstitute.org";

                log.error(msg);
                MessageUtils.showMessage(msg);
                break;
            }
            lastAllocatedCount = allocatedCount;

            // Loop through alignments until we reach the end of the interval

            PriorityQueue<T> bucket = null;

            ArrayList<Integer> emptyBucketKeys = new ArrayList();
            for (Integer key : bucketArray.keySet()) {
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
                lastAllocatedCount = -1;
            }
            currentRow = new FeatureRow();
            nextStart = 0;
        }
        // Add the last row
        if (currentRow.features.size() > 0) {
            rows.add(currentRow);
        }

        return rows;
    }

    /**
     * Takes in an iterator of Features and creates two new iterators, one that gives the neg-strand features
     * and one that gives the pos-strand features.  Does this by buffering features of the opposite strand
     * when a feature of a particular strand is requested.
     */
    protected class IteratorSplitterByCharge {
        Iterator origIter;

        List<Feature> posBuffer = new ArrayList<Feature>();
        List<Feature> negBuffer = new ArrayList<Feature>();

        PerStrandIter posIter;
        PerStrandIter negIter;

        public IteratorSplitterByCharge(Iterator origIter) {
            this.origIter = origIter;

            posIter = new PerStrandIter(posBuffer, negBuffer, Strand.POSITIVE);
            negIter = new PerStrandIter(negBuffer, posBuffer, Strand.NEGATIVE);
        }

        public PerStrandIter getPosIter() {
            return posIter;
        }

        public PerStrandIter getNegIter() {
            return negIter;
        }

        /**
         * An iterator of Features for a single strand. Polls the original iterator for new features and
         * either returns them or buffers them in the other strand's buffer, accordingly
         */
        protected class PerStrandIter implements Iterator<Feature> {
            List<Feature> buffer;
            List<Feature> otherBuffer;
            Strand strand;


            public PerStrandIter(List<Feature> buffer, List<Feature> otherBuffer, Strand strand) {
                this.buffer = buffer;
                this.otherBuffer = otherBuffer;
                this.strand = strand;
            }

            /**
             * Burn through the original iterator until we either hit the end or find a feature of the
             * appropriate strand, buffering the other strand features as we go
             * @return
             */
            public boolean hasNext() {
                while (buffer.isEmpty() && origIter.hasNext())
                {
                    BasicFeature feature = (BasicFeature) origIter.next();
                    if (feature.getStrand() == strand)
                        buffer.add(feature);
                    else
                        otherBuffer.add(feature);
                }
                return !buffer.isEmpty();
            }

            /**
             * Burn through the original iterator until we either hit the end or find a feature of the
             * appropriate strand, buffering the other strand features as we go
             * @return
             */
            public Feature next() {
                while (buffer.isEmpty() && origIter.hasNext())
                {
                    BasicFeature feature = (BasicFeature) origIter.next();
                    if (feature.getStrand() == strand)
                        buffer.add(feature);
                    else
                        otherBuffer.add(feature);
                }
                if (buffer.isEmpty())
                    return null;
                Feature result = buffer.get(0);
                buffer.remove(0);
                return result;
            }

            public void remove() {
                //not implemented
            }
        }
    }

}
