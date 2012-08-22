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

import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.PreferenceManager;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.sam.reader.ReadGroupFilter;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Ignore;
import org.junit.Test;

import java.io.IOException;
import java.util.*;

import static org.junit.Assert.*;

/**
 * @author jrobinso
 */
public class AlignmentIntervalLoaderTest extends AbstractHeadlessTest {

    String testFile = "http://www.broadinstitute.org/igvdata/BodyMap/hg18/50bp/FCA/s_1_1_sequence.bam";
    String sequence = "chr1";
    int start = 44680145;
    int end = 44789983;
    private boolean contained = false;


    @Test
    public void testQuery() throws IOException {
        tstQuery(testFile, sequence, start, end, contained, Integer.MAX_VALUE / 1000);
    }

    /**
     * Test of query method, of class AlignmentIntervalLoader.  The test compares
     * the results of AlignmentIntervalLoader  with an AlignmentReader.
     * <p/>
     * Note that SAMFileReader (which is the non-caching reader) is 1-based
     * and inclusive-end. AlignmentIntervalLoader is 0-based and exclusive end.
     */
    public void tstQuery(String testFile, String sequence, int start, int end, boolean contained, int maxDepth) throws IOException {

        ResourceLocator loc = new ResourceLocator(testFile);
        AlignmentReader reader = AlignmentReaderFactory.getReader(loc);
        CloseableIterator<Alignment> iter = reader.query(sequence, start, end, contained);

        List<Alignment> expectedResult = new ArrayList<Alignment>();
        while (iter.hasNext()) {
            Alignment rec = iter.next();

            // the following filters are applied in the Caching reader, so we need to apply them here.
            boolean filterFailedReads = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_FILTER_FAILED_READS);
            ReadGroupFilter filter = ReadGroupFilter.getFilter();
            boolean showDuplicates = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_SHOW_DUPLICATES);
            int qualityThreshold = PreferenceManager.getInstance().getAsInt(PreferenceManager.SAM_QUALITY_THRESHOLD);
            if (!rec.isMapped() || (!showDuplicates && rec.isDuplicate()) ||
                    (filterFailedReads && rec.isVendorFailedRead()) ||
                    rec.getMappingQuality() < qualityThreshold ||
                    (filter != null && filter.filterAlignment(rec))) {
                continue;
            }

            expectedResult.add(rec);
            //System.out.println("name: " + rec.getReadName() + "strt: " + rec.getStart() + " end: " + rec.getEnd());
            if (contained) {
                assertTrue(rec.getStart() >= start);
            } else {
                //All we require is some overlap
                boolean overlap = rec.getStart() >= start && rec.getStart() < end;
                overlap |= (rec.getEnd() >= start) && (rec.getStart() < start);
                assertTrue(overlap);
            }
            assertEquals(sequence, rec.getChr());
        }
        reader.close();

        reader = AlignmentReaderFactory.getReader(loc);
        AlignmentIntervalLoader alignmentIntervalLoader = new AlignmentIntervalLoader(reader);

        boolean showSpliceJunctions = false;
        AlignmentTrack.RenderOptions renderOptions = new AlignmentTrack.RenderOptions();
        Map<String, PEStats> peStats = new HashMap<String, PEStats>();
        AlignmentTrack.BisulfiteContext bisulfiteContext = null;

        AlignmentInterval interval = alignmentIntervalLoader.loadInterval(sequence, start, end, showSpliceJunctions,
                renderOptions, peStats, bisulfiteContext);
        List<Alignment> result = new ArrayList();

        Iterator<Alignment> alignmentIterator = interval.getAlignmentIterator();
        while (alignmentIterator.hasNext()) {
            result.add(alignmentIterator.next());
        }


        assertTrue(expectedResult.size() > 0);
        assertEquals(expectedResult.size(), result.size());

        //Reads sorted by start position, apparently there is some wiggle room in the exact order
        //We sort each first by start position, then end position
        Collections.sort(expectedResult, new StartEndSorter());
        Collections.sort(result, new StartEndSorter());
        for (int i = 0; i < result.size(); i++) {
            Alignment rec = result.get(i);

            if (i % 2 == 0 && rec.isPaired()) {
                //Check that paired reads are together

                System.out.println(rec.getReadName());

                System.out.println(result.get(i + 1).getReadName());
                //assertEquals(rec.getReadName(), result.get(i+1).getReadName());
            }

            if (contained) {
                assertTrue(rec.getStart() >= start);
            } else {
                //All we require is some overlap
                boolean overlap = rec.getStart() >= start && rec.getStart() < end;
                overlap |= start >= rec.getStart() && start < rec.getEnd();
                assertTrue(overlap);
            }
            assertEquals(sequence, rec.getChr());

            Alignment exp = expectedResult.get(i);
            assertEquals("Start mismatch at position " + i + " read name " + exp.getReadName(), exp.getStart(), rec.getStart());

            assertEquals(exp.getReadName(), rec.getReadName());
            assertEquals("End mismatch at position " + i + " read name " + rec.getReadName(), exp.getEnd(), rec.getEnd());
        }
    }

    @Ignore
    @Test
    public void testQueryLargeFile() throws Exception {
        PreferenceManager.getInstance().put(PreferenceManager.SAM_MAX_VISIBLE_RANGE, "5");
        String path = TestUtils.LARGE_DATA_DIR + "ABCD_igvSample.bam";

        ResourceLocator loc = new ResourceLocator(path);
        AlignmentReader reader = AlignmentReaderFactory.getReader(loc);
        AlignmentIntervalLoader alignmentIntervalLoader = new AlignmentIntervalLoader(reader);

        //Edge location
        String sequence = "chr12";
        int start = 56815621 - 1;
        int end = start + 1;
        int expSize = 1066;

        tstSize(alignmentIntervalLoader, sequence, start, end, (int) (expSize * 1.6), expSize);
        tstQuery(path, sequence, start, end, false, 10000);

        //Edge location, downsampled
        sequence = "chr12";
        start = 56815635 - 1;
        end = start + 1;
        expSize = 165;

        reader = AlignmentReaderFactory.getReader(loc);
        alignmentIntervalLoader = new AlignmentIntervalLoader(reader);

        tstSize(alignmentIntervalLoader, sequence, start, end, expSize + 20, expSize);
        tstQuery(path, sequence, start, end, false, 10000);

        //Center location
        sequence = "chr12";
        start = 56815675 - 1;
        end = start + 1;

        expSize = 3288;

        reader = AlignmentReaderFactory.getReader(loc);
        alignmentIntervalLoader = new AlignmentIntervalLoader(reader);

        tstSize(alignmentIntervalLoader, sequence, start, end, expSize + 20, expSize);
        tstQuery(path, sequence, start, end, false, 10000);
    }

    @Test
    public void testQueryPiledUp() throws Exception {
        PreferenceManager.getInstance().put(PreferenceManager.SAM_MAX_VISIBLE_RANGE, "5");
        PreferenceManager.getInstance().put(PreferenceManager.SAM_DOWNSAMPLE_READS, "false");
        String path = TestUtils.DATA_DIR + "aligned/pileup.sorted.aligned";

        ResourceLocator loc = new ResourceLocator(path);
        AlignmentReader reader = AlignmentReaderFactory.getReader(loc);
        AlignmentIntervalLoader alignmentIntervalLoader = new AlignmentIntervalLoader(reader);

        //Edge location
        String sequence = "chr1";
        int start = 141 - 1;
        int end = start + 1;
        int expSize = 40;

        tstSize(alignmentIntervalLoader, sequence, start, end, expSize * 7, expSize);

        loc = new ResourceLocator(path);
        reader = AlignmentReaderFactory.getReader(loc);
        alignmentIntervalLoader = new AlignmentIntervalLoader(reader);

        tstSize(alignmentIntervalLoader, sequence, start, end, expSize * 100, expSize);
        tstQuery(path, sequence, start, end, false, 10000);

        //Center, deep coverage region
        sequence = "chr1";
        start = 429;
        end = start + 1;
        int coverageLim = 1000;
        expSize = 1408;

        loc = new ResourceLocator(path);
        reader = AlignmentReaderFactory.getReader(loc);
        alignmentIntervalLoader = new AlignmentIntervalLoader(reader);


        tstSize(alignmentIntervalLoader, sequence, start, end, coverageLim, expSize);

        coverageLim = 10000;
        expSize = 1408;

        loc = new ResourceLocator(path);
        reader = AlignmentReaderFactory.getReader(loc);
        alignmentIntervalLoader = new AlignmentIntervalLoader(reader);


        tstSize(alignmentIntervalLoader, sequence, start, end, coverageLim, expSize);

        tstQuery(path, sequence, start, end, false, coverageLim);

    }

    public List<Alignment> tstSize(AlignmentIntervalLoader alignmentIntervalLoader, String sequence, int start, int end, int maxDepth, int expSize) {

        boolean showSpliceJunctions = false;
        AlignmentTrack.RenderOptions renderOptions = new AlignmentTrack.RenderOptions();
        AlignmentDataManager.DownsampleOptions downsampleOptions = new AlignmentDataManager.DownsampleOptions();
        Map<String, PEStats> peStats = new HashMap<String, PEStats>();
        AlignmentTrack.BisulfiteContext bisulfiteContext = null;

        AlignmentInterval interval = alignmentIntervalLoader.loadInterval(sequence, start, end, showSpliceJunctions,
                renderOptions, peStats, bisulfiteContext);
        List<Alignment> result = new ArrayList();

        Iterator<Alignment> iter = interval.getAlignmentIterator();
        while (iter.hasNext()) {
            result.add(iter.next());
        }

        assertEquals(expSize, result.size());
        return result;
    }


    /**
     * The main purpose of this test is to see if we get a
     * heap space error.
     *
     * @throws Exception
     */
    @Test
    public void testQueryLargeFile2() throws Exception {
        String path = "http://www.broadinstitute.org/igvdata/1KG/pilot2Bams/NA12878.454.bam";

        ResourceLocator loc = new ResourceLocator(path);
        AlignmentReader reader = AlignmentReaderFactory.getReader(loc);
        AlignmentIntervalLoader alignmentIntervalLoader = new AlignmentIntervalLoader(reader);

        String sequence = "MT";
        int start = 1000;
        int end = 3000;


        boolean showSpliceJunctions = false;
        AlignmentTrack.RenderOptions renderOptions = new AlignmentTrack.RenderOptions();
        Map<String, PEStats> peStats = new HashMap<String, PEStats>();
        AlignmentTrack.BisulfiteContext bisulfiteContext = null;

        AlignmentInterval interval = alignmentIntervalLoader.loadInterval(sequence, start, end, showSpliceJunctions,
                renderOptions, peStats, bisulfiteContext);

        Iterator<Alignment> iter = interval.getAlignmentIterator();
        int count = 0;
        while (iter.hasNext()) {
            count++;
        }

        assertTrue(count > 0);

    }

    /**
     * Sorts by: read name, start, end
     */
    private class StartEndSorter implements Comparator<Alignment> {

        public int compare(Alignment o1, Alignment o2) {
            Alignment al1 = (Alignment) o1;
            Alignment al2 = (Alignment) o2;

            String n1 = al1.getReadName();
            n1 = n1 != null ? n1 : "";
            int nStart = n1.compareTo(al2.getReadName());
            if (nStart != 0) {
                return nStart;
            }

            int cStarts = compareInts(al1.getStart(), al2.getStart());
            if (cStarts != 0) {
                return cStarts;
            }

            int cEnds = compareInts(al1.getEnd(), al2.getEnd());
            return cEnds;
        }

        private int compareInts(int i1, int i2) {
            if (i1 < i2) {
                return -1;
            } else if (i1 > i2) {
                return +1;
            } else {
                return 0;
            }

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