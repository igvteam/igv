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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.PreferenceManager;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.util.ResourceLocator;
import org.junit.Ignore;
import org.junit.Test;

import java.util.*;

import static org.junit.Assert.*;

/**
 * @author jrobinso
 */
public class AlignmentTileLoaderTest extends AbstractHeadlessTest {


    /**
     * Test that sampling keeps pairs together. Note that this test requires a lot of memory (2Gb on this devs machine)
     * Pairs are only definitely kept together if they end up on the same tile, so we need 1 big tile.
     *
     * @throws Exception
     */
    @Ignore
    @Test
    public void testKeepPairs() throws Exception {
        String path = "http://www.broadinstitute.org/igvdata/1KG/pilot2Bams/NA12878.SOLID.bam";


        String sequence = "1";
        int start = 1;
        int end = 2000;
        int maxDepth = 2;
        String max_vis = PreferenceManager.getInstance().get(PreferenceManager.SAM_MAX_VISIBLE_RANGE);
        PreferenceManager.getInstance().put(PreferenceManager.SAM_MAX_VISIBLE_RANGE, "" + (end - start));

        try {
            ResourceLocator loc = new ResourceLocator(path);
            AlignmentReader reader = AlignmentReaderFactory.getReader(loc);
            AlignmentTileLoader loader = new AlignmentTileLoader(reader);

            AlignmentDataManager.DownsampleOptions downsampleOptions = new AlignmentDataManager.DownsampleOptions();

            AlignmentTileLoader.AlignmentTile tile = loader.loadTile(sequence, start, end, null, downsampleOptions, null, null, null);
            List<Alignment> alignments = tile.getAlignments();
            int count = 0;
            Map<String, Integer> pairedReads = new HashMap<String, Integer>();
            for(Alignment al: alignments) {
                assertNotNull(al);
                count++;

                if (al.isPaired() && al.getMate().isMapped()) {
                    //Mate may not be part of the query.
                    //Make sure it's within bounds
                    int mateStart = al.getMate().getStart();
                    //All we require is some overlap
                    boolean overlap = mateStart >= start && mateStart < end;
                    if (overlap) {
                        Integer rdCnt = pairedReads.get(al.getReadName());
                        rdCnt = rdCnt != null ? rdCnt + 1 : 1;
                        pairedReads.put(al.getReadName(), rdCnt);
                    }
                }
            }

            assertTrue(count > 0);

            //Note: CachingQueryReader will not line alignments up properly if they land in different tiles
            int countmissing = 0;
            for (String readName : pairedReads.keySet()) {
                int val = pairedReads.get(readName);
                countmissing += 2 == val ? 0 : 1;
                if (val != 2) {
                    System.out.println("Read " + readName + " has val " + val);
                }
            }

            //System.out.println("Number of paired reads: " + pairedReads.size());
            assertTrue("No pairs in test data set", pairedReads.size() > 0);
            assertEquals("Missing " + countmissing + " out of " + pairedReads.size() + " pairs", 0, countmissing);
        } catch (Exception e) {
            PreferenceManager.getInstance().put(PreferenceManager.SAM_MAX_VISIBLE_RANGE, max_vis);
            throw e;
        }

    }


    public static void main(String[] args) {
        //Represents total number of alignments
        long totalLength = (long) 1e6;
        //Memory used per alignment
        int longseach = 100;

        int maxKeep = 1000;
        float fmaxKeep = (float) maxKeep;
        int maxBucketDepth = (int) 1e5;

        long seed = 5310431327l;

        long t1 = System.currentTimeMillis();
        liveSample(totalLength, longseach, seed, maxKeep, maxBucketDepth * 10);
        long t2 = System.currentTimeMillis();
        System.out.println("Time for live sampling: " + (t2 - t1) + " mSec");

        long t3 = System.currentTimeMillis();
        downSample(totalLength, longseach, seed, maxKeep, maxBucketDepth);
        long t4 = System.currentTimeMillis();
        System.out.println("Time for down sampling: " + (t4 - t3) + " mSec");

    }

    // NOTE: The active tests in this class have been moved to AlignmentDataManagerTest.

    @Test
    public void nullTest() {
        assert(true == true);
    }

    /**
     * Test that our live sample gives a uniform distribution
     */
    @Ignore
    public void testLiveSample() throws Exception {
        int totalLength = (int) 1e4;
        //Store the number of times each index is sampled
        int[] counts = new int[totalLength];
        List<long[]> samples;
        int longseach = 1;
        long seed = 212338399;
        Random rand = new Random(seed);
        int maxKeep = 1000;
        int maxBucketDepth = Integer.MAX_VALUE;

        int trials = 10000;
        for (int _ = 0; _ < trials; _++) {
            seed = rand.nextLong();
            samples = liveSample(totalLength, longseach, seed, maxKeep, maxBucketDepth);
            for (long[] dat : samples) {
                counts[(int) dat[0]] += 1;
            }
        }

        float avgFreq = ((float) maxKeep) / totalLength;
        int avgCount = (int) (avgFreq * trials);
        double stdDev = Math.sqrt(trials / 12);
        int numStds = 4;

        int ind = 0;
        //System.out.println("Expected number of times sampled: " + avgCount + ". Stdev " + stdDev);
        for (int cnt : counts) {
            //System.out.println("ind: " + ind + " cnt: " + cnt);
            assertTrue("Index " + ind + " outside of expected sampling range at " + cnt, Math.abs(cnt - avgCount) < numStds * stdDev);
            ind++;
        }


    }

    private static List<long[]> liveSample(long totalLength, int longseach, long seed, int maxKeep, int maxBucketDepth) {

        List<long[]> liveSampled = new ArrayList<long[]>(maxKeep);
        float fmaxKeep = (float) maxKeep;

        RandDataIterator iter1 = new RandDataIterator(totalLength, longseach);
        float prob = 0;
        Random rand = new Random(seed);
        int numAfterMax = 1;
        for (long[] data : iter1) {
            if (liveSampled.size() < maxKeep) {
                liveSampled.add(data);
            } else if (liveSampled.size() > maxBucketDepth) {
                break;
            } else {
                //Calculate whether to accept this element
                prob = fmaxKeep / (maxKeep + numAfterMax);
                numAfterMax += 1;
                boolean keep = rand.nextFloat() < prob;
                if (keep) {
                    //Choose one to replace
                    int torep = rand.nextInt(maxKeep);
                    liveSampled.remove(torep);
                    liveSampled.add(data);
                }

            }
        }
        return liveSampled;
    }

    private static List<long[]> downSample(long totalLength, int longseach, long seed, int maxKeep, int maxBucketDepth) {

        List<long[]> downSampled = new ArrayList<long[]>(maxKeep);
        RandDataIterator iter2 = new RandDataIterator(totalLength, longseach);
        Random rand = new Random(seed);
        for (long[] data : iter2) {
            if (downSampled.size() < maxBucketDepth) {
                downSampled.add(data);
            } else {
                break;
            }
        }

        //Actual downsampling
        while (downSampled.size() > maxKeep) {
            downSampled.remove(rand.nextInt(downSampled.size()));
        }

        return downSampled;
    }

    /**
     * Iterator over garbage data.
     */
    private static class RandDataIterator implements Iterable<long[]>, Iterator<long[]> {

        private long counter;
        private long length;
        private int longseach;

        /**
         * @param length    Number of elements that this iterator will have
         * @param longseach how large each element will be (byte array)
         */
        public RandDataIterator(long length, int longseach) {
            this.length = length;
            this.longseach = longseach;
        }

        public boolean hasNext() {
            return counter < length;
        }

        public long[] next() {
            if (!hasNext()) {
                return null;
            }
            long[] arr = new long[longseach];
            Arrays.fill(arr, counter);
            counter++;
            return arr;
        }

        public void remove() {
            throw new UnsupportedOperationException("Can't remove");
        }

        public Iterator<long[]> iterator() {
            return this;
        }
    }


}