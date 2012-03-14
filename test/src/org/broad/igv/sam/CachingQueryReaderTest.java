/*
 * Copyright (c) 2007-2012 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTIES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NON-INFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
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

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.PreferenceManager;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.sam.reader.ReadGroupFilter;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.*;

import java.io.IOException;
import java.util.*;

import static org.junit.Assert.*;

/**
 * @author jrobinso
 */
public class CachingQueryReaderTest {

    String testFile = "http://www.broadinstitute.org/igvdata/BodyMap/hg18/50bp/FCA/s_1_1_sequence.bam";
    String sequence = "chr1";
    int start = 44680145;
    int end = 44789983;
    private boolean contained = false;

    public CachingQueryReaderTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        TestUtils.setUpHeadless();
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of getHeader method, of class CachingQueryReader. The test compares
     * the results of CachingQueryReader with a non-caching reader which
     * is assumed to be correct.
     */
    @Test
    public void testGetHeader() throws IOException {

        ResourceLocator loc = new ResourceLocator(testFile);
        AlignmentReader reader = AlignmentReaderFactory.getReader(loc);
        SAMFileHeader expectedHeader = reader.getHeader();
        reader.close();

        reader = AlignmentReaderFactory.getReader(loc);
        CachingQueryReader cachingReader = new CachingQueryReader(reader);
        SAMFileHeader header = cachingReader.getHeader();
        cachingReader.close();

        assertTrue(header.equals(expectedHeader));
    }

    @Test
    public void testQuery() throws IOException {
        tstQuery(testFile, sequence, start, end, contained, Integer.MAX_VALUE / 1000);
    }

    /**
     * Test of query method, of class CachingQueryReader.  The test compares
     * the results of CachingQueryReader non-caching reader which
     * is assumed to be correct.
     * <p/>
     * Note that SAMFileReader (which is the non-caching reader) is 1-based
     * and inclusive-end. CachingQueryReader is 0-based and exclusive end.
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
        CachingQueryReader cachingReader = new CachingQueryReader(reader);
        CloseableIterator<Alignment> cachingIter = cachingReader.query(sequence, start, end, new ArrayList(),
                new ArrayList(), maxDepth, null, null);
        List<Alignment> result = new ArrayList();

        while (cachingIter.hasNext()) {
            result.add(cachingIter.next());
        }
        cachingReader.close();


        assertTrue(expectedResult.size() > 0);
        assertEquals(expectedResult.size(), result.size());

        //Reads sorted by start position, apparently there is some wiggle room in the exact order
        //We sort each first by start position, then end position
        Collections.sort(expectedResult, new StartEndSorter());
        Collections.sort(result, new StartEndSorter());
        for (int i = 0; i < result.size(); i++) {
            Alignment rec = result.get(i);

            if (contained) {
                assertTrue(rec.getStart() >= start);
            } else {
                //All we require is some overlap
                boolean overlap = rec.getStart() >= start && rec.getStart() <= end;
                overlap |= start >= rec.getStart() && start <= rec.getEnd();
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
        String path = TestUtils.LARGE_DATA_DIR + "/ABCD_igvSample.bam";

        ResourceLocator loc = new ResourceLocator(path);
        AlignmentReader reader = AlignmentReaderFactory.getReader(loc);
        CachingQueryReader cachingReader = new CachingQueryReader(reader);

        //Edge location
        String sequence = "chr12";
        int start = 56815621 - 1;
        int end = start + 1;
        int expSize = 1066;

        tstSize(cachingReader, sequence, start, end, expSize * 5, expSize);
        tstQuery(path, sequence, start, end, false, 10000);

        //Edge location, downsampled
        sequence = "chr12";
        start = 56815635 - 1;
        end = start + 1;
        expSize = 165;

        tstSize(cachingReader, sequence, start, end, expSize + 2, expSize);
        tstQuery(path, sequence, start, end, false, 10000);

        //Center location
        sequence = "chr12";
        start = 56815674;
        end = start + 1;

        expSize = 3288;

        tstSize(cachingReader, sequence, start, end, expSize * 5, expSize);
        tstQuery(path, sequence, start, end, false, 10000);


    }

    @Test
    public void testQueryPiledUp() throws Exception {
        PreferenceManager.getInstance().put(PreferenceManager.SAM_MAX_VISIBLE_RANGE, "5");
        String path = TestUtils.DATA_DIR + "/aligned/pileup.sorted.aligned";

        ResourceLocator loc = new ResourceLocator(path);
        AlignmentReader reader = AlignmentReaderFactory.getReader(loc);
        CachingQueryReader cachingReader = new CachingQueryReader(reader);

        //Edge location
        String sequence = "chr1";
        int start = 141 - 1;
        int end = start + 1;
        int expSize = 40;

        tstSize(cachingReader, sequence, start, end, expSize * 7, expSize);

        loc = new ResourceLocator(path);
        reader = AlignmentReaderFactory.getReader(loc);
        cachingReader = new CachingQueryReader(reader);

        tstSize(cachingReader, sequence, start, end, expSize * 100, expSize);
        //tstQuery(path, sequence, start, end, false, 10000);

        //Center, deep coverage region
        sequence = "chr1";
        start = 429;
        end = start + 1;
        int coverageLim = 1000;
        expSize = 1408;

        loc = new ResourceLocator(path);
        reader = AlignmentReaderFactory.getReader(loc);
        cachingReader = new CachingQueryReader(reader);


        tstSize(cachingReader, sequence, start, end, coverageLim, expSize);

        coverageLim = 10000;
        expSize = 1408;

        loc = new ResourceLocator(path);
        reader = AlignmentReaderFactory.getReader(loc);
        cachingReader = new CachingQueryReader(reader);


        tstSize(cachingReader, sequence, start, end, coverageLim, expSize);

        tstQuery(path, sequence, start, end, false, coverageLim);

    }

    public List<Alignment> tstSize(CachingQueryReader cachingReader, String sequence, int start, int end, int maxDepth, int expSize) {
        CloseableIterator<Alignment> cachingIter = cachingReader.query(sequence, start, end, new ArrayList(),
                new ArrayList(), maxDepth, null, null);
        List<Alignment> result = new ArrayList();

        while (cachingIter.hasNext()) {
            result.add(cachingIter.next());
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
        CachingQueryReader cachingReader = new CachingQueryReader(reader);

        String sequence = "MT";
        int start = 1000;
        int end = 3000;
        int maxDepth = 500;

        CloseableIterator<Alignment> iter = cachingReader.query(sequence, start, end, new ArrayList(),
                new ArrayList(), maxDepth, null, null);
        int count = 0;
        while (iter.hasNext()) {
            assertNotNull(iter.next());
            count++;
        }

        assertTrue(count > 0);


    }

    private class StartEndSorter implements Comparator<Alignment> {

        @Override
        public int compare(Alignment o1, Alignment o2) {
            Alignment al1 = (Alignment) o1;
            Alignment al2 = (Alignment) o2;

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
    @Test
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

        @Override
        public boolean hasNext() {
            return counter < length;
        }

        @Override
        public long[] next() {
            if (!hasNext()) {
                return null;
            }
            long[] arr = new long[longseach];
            Arrays.fill(arr, counter);
            counter++;
            return arr;
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("Can't remove");
        }

        @Override
        public Iterator<long[]> iterator() {
            return this;
        }
    }


}