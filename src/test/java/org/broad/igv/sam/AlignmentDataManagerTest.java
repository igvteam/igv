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

import htsjdk.samtools.util.CloseableIterator;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.sam.reader.ReadGroupFilter;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

import java.io.IOException;
import java.util.*;

import static junit.framework.Assert.assertNotNull;
import static junit.framework.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012-Jul-12
 */
public class AlignmentDataManagerTest extends AbstractHeadlessTest {

    @Rule
    public TestRule testTimeout = new Timeout((int) 120e3);

    private static String frameName = "testFrame";

    @Override
    public void tearDown() throws Exception {
        super.tearDown();
    }
//
//
//    @Test
//    public void testPreloadNoMerge() throws Exception {
//
//        final String chr = "chr1";
//        final int start = 151666494;
//        final int halfwidth = 1000;
//        final int end = start + 2 * halfwidth;
//
//        //Load separate intervals, check they don't merge
//        AlignmentDataManager manager = getManager171();
//        ReferenceFrame frame = new ReferenceFrame(frameName);
//        AlignmentTrack.RenderOptions renderOptions = new AlignmentTrack.RenderOptions();
//        frame.setBounds(0, end - start);
//
//        RenderContext context = new RenderContext(null, null, frame, null);
//
//        int lastStart = genome.getChromosome(chr).getLength() - 4 * halfwidth;
//        int[] starts = new int[]{500, 5000, 15000, start, 500000, lastStart};
//        int[] ends = new int[]{600, 10000, 20000, end, 600000, lastStart + 2 * halfwidth};
//        for (int ii = 0; ii < starts.length; ii++) {
//            frame.jumpTo(new Locus(chr, starts[ii], ends[ii]));
//            int actEnd = (int) frame.getEnd();
//
//            manager.load(context.getReferenceFrame(), renderOptions, false);
//
//            assertManagerHasInterval(manager, context.getReferenceFrame(), chr, starts[ii], actEnd);
//        }
//
//
//    }

    private static void assertManagerHasInterval(AlignmentDataManager manager, ReferenceFrame frame, String chr, int start, int end) {

        AlignmentInterval interval = manager.getLoadedInterval(frame);
        assertNotNull(interval);

        boolean haveInterval = interval.contains(chr, start, end);
        assertTrue(haveInterval);
    }

    public static AlignmentDataManager getManager171() throws IOException {

        String infilepath = TestUtils.LARGE_DATA_DIR + "HG00171.hg18.bam";
        ResourceLocator locator = new ResourceLocator(infilepath);
        AlignmentDataManager manager = new AlignmentDataManager(locator, genome);
        return manager;
    }



    @Test
    public void testQuery() throws IOException {
        String testFile = TestUtils.DATA_DIR + "bam/gstt1_sample.bam";
        String sequence = "chr22";
        int start = 24376039;
        int end = 24376625;
        boolean contained = false;

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
            boolean filterFailedReads = PreferencesManager.getPreferences().getAsBoolean(Constants.SAM_FILTER_FAILED_READS);
            ReadGroupFilter filter = ReadGroupFilter.getFilter();
            boolean showDuplicates = !PreferencesManager.getPreferences().getAsBoolean(Constants.SAM_FILTER_DUPLICATES);
            int qualityThreshold = PreferencesManager.getPreferences().getAsInt(Constants.SAM_QUALITY_THRESHOLD);
            if (!rec.isMapped() || (!showDuplicates && rec.isDuplicate()) ||
                    (filterFailedReads && rec.isVendorFailedRead()) ||
                    rec.getMappingQuality() < qualityThreshold ||
                    (filter != null && filter.filterAlignment(rec))) {
                continue;
            }

            expectedResult.add(rec);
            //System.out.println("name: " + rec.getReadName() + "strt: " + rec.getStart() + " end: " + rec.getEnd());
            if (contained) {
                Assert.assertTrue(rec.getStart() >= start);
            } else {
                //All we require is some overlap
                boolean overlap = rec.getStart() >= start && rec.getStart() < end;
                overlap |= (rec.getEnd() >= start) && (rec.getStart() < start);
                Assert.assertTrue(overlap);
            }
            Assert.assertEquals(sequence, rec.getChr());
        }
        reader.close();
        AlignmentDataManager manager = new AlignmentDataManager(loc, genome);

        // Turn off downsampling
        IGVPreferences prefs = PreferencesManager.getPreferences();
        prefs.put(Constants.SAM_DOWNSAMPLE_READS, "false");

        AlignmentInterval interval = loadInterval(manager, sequence, start, end);
        List<Alignment> result = new ArrayList();

        Iterator<Alignment> alignmentIterator = interval.getAlignmentIterator();
        while (alignmentIterator.hasNext()) {
            result.add(alignmentIterator.next());
        }


        Assert.assertTrue(expectedResult.size() > 0);
        Assert.assertEquals(expectedResult.size(), result.size());

        //Reads sorted by start position, apparently there is some wiggle room in the exact order
        //We sort each first by start position, then end position
        Collections.sort(expectedResult, new StartEndSorter());
        Collections.sort(result, new StartEndSorter());
        for (int i = 0; i < result.size(); i++) {
            Alignment rec = result.get(i);
            
            if (contained) {
                Assert.assertTrue(rec.getStart() >= start);
            } else {
                //All we require is some overlap
                boolean overlap = rec.getStart() >= start && rec.getStart() < end;
                overlap |= start >= rec.getStart() && start < rec.getEnd();
                Assert.assertTrue(overlap);
            }
            Assert.assertEquals(sequence, rec.getChr());

            Alignment exp = expectedResult.get(i);
            Assert.assertEquals("Start mismatch at position " + i + " read name " + exp.getReadName(), exp.getStart(), rec.getStart());

            Assert.assertEquals(exp.getReadName(), rec.getReadName());
            Assert.assertEquals("End mismatch at position " + i + " read name " + rec.getReadName(), exp.getEnd(), rec.getEnd());
        }
    }

    @Ignore
    @Test
    public void testQueryLargeFile() throws Exception {
        PreferencesManager.getPreferences().put(Constants.SAM_MAX_VISIBLE_RANGE, "5");
        String path = TestUtils.LARGE_DATA_DIR + "ABCD_igvSample.bam";

        ResourceLocator loc = new ResourceLocator(path);
        AlignmentDataManager manager = new AlignmentDataManager(loc, genome);

        //Edge location
        String sequence = "chr12";
        int start = 56815621 - 1;
        int end = start + 1;
        int expSize = 1066;

        tstSize(manager, sequence, start, end, (int) (expSize * 1.6), expSize);
        tstQuery(path, sequence, start, end, false, 10000);

        //Edge location, downsampled
        sequence = "chr12";
        start = 56815635 - 1;
        end = start + 1;
        expSize = 165;


        tstSize(manager, sequence, start, end, expSize + 20, expSize);
        tstQuery(path, sequence, start, end, false, 10000);

        //Center location
        sequence = "chr12";
        start = 56815675 - 1;
        end = start + 1;

        expSize = 3288;

        tstSize(manager, sequence, start, end, expSize + 20, expSize);
        tstQuery(path, sequence, start, end, false, 10000);
    }

    @Test @Ignore("Fails unless tests are run in separate JVMs")
    public void testQueryPiledUp() throws Exception {
        PreferencesManager.getPreferences().put(Constants.SAM_MAX_VISIBLE_RANGE, "5");
        PreferencesManager.getPreferences().put(Constants.SAM_DOWNSAMPLE_READS, "false");
        String path = TestUtils.DATA_DIR + "aligned/pileup.sorted.aligned";

        ResourceLocator loc = new ResourceLocator(path);
        AlignmentDataManager manager = new AlignmentDataManager(loc, genome);

        //Edge location
        String sequence = "chr1";
        int start = 141 - 1;
        int end = start + 1;
        int expSize = 40;

        tstSize(manager, sequence, start, end, expSize * 7, expSize);

        tstSize(manager, sequence, start, end, expSize * 100, expSize);
        tstQuery(path, sequence, start, end, false, 10000);

        //Center, deep coverage region
        sequence = "chr1";
        start = 429;
        end = start + 1;
        int coverageLim = 1000;
        expSize = 1408;

        tstSize(manager, sequence, start, end, coverageLim, expSize);

        coverageLim = 10000;
        expSize = 1408;

        tstSize(manager, sequence, start, end, coverageLim, expSize);

        tstQuery(path, sequence, start, end, false, coverageLim);

    }

    public List<Alignment> tstSize(AlignmentDataManager manager, String sequence, int start, int end, int maxDepth, int expSize) {
        AlignmentInterval interval = loadInterval(manager, sequence, start, end);

        List<Alignment> result = new ArrayList();

        Iterator<Alignment> iter = interval.getAlignmentIterator();
        while (iter.hasNext()) {
            result.add(iter.next());
        }

        Assert.assertEquals(expSize, result.size());
        return result;
    }


    /**
     * The main purpose of this test is to see if we get a
     * heap space error.
     *
     * @throws Exception
     */
    @Ignore
    @Test
    public void testQueryLargeFile2() throws Exception {
        String path = "http://1000genomes.s3.amazonaws.com/data/NA12878/high_coverage_alignment/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130520.bam";

        ResourceLocator loc = new ResourceLocator(path);

        String sequence = "chrM";
        int start = 0;
        int end = 200;

        System.gc();

        long startTime = System.nanoTime();
        AlignmentDataManager manager = new AlignmentDataManager(loc, genome);

        AlignmentInterval interval = loadInterval(manager, sequence, start, end);

        System.out.println("# of downsampled intervals: " + interval.getDownsampledIntervals().size());
        Iterator<Alignment> iter = interval.getAlignmentIterator();

        int count = 0;
        while (iter.hasNext()) {
            count++;
            Alignment al = iter.next();
            assertTrue(al.getStart() <= end);
            assertTrue(al.getEnd() >= start);
        }
        long endTime = System.nanoTime();
        long total = endTime - startTime;
        System.out.println(String.format("Total time: %2.2f sec, %d alignments kept\"", total * 1.0 / 1e9, count));

        Assert.assertTrue(count > 0);
    }

    /**
     * Very basic test that nothing crashes when loading an alignment
     * which is padded
     * @throws Exception
     */
    @Test
    public void testWithPadding() throws Exception{
        String filepath = TestUtils.DATA_DIR + "sam/has_padding.sam";
        TestUtils.createIndex(filepath);

        AlignmentDataManager manager = new AlignmentDataManager(new ResourceLocator(filepath), genome);
        AlignmentInterval interval = manager.loadInterval("chr22", 0, Integer.MAX_VALUE, null);
        Iterator<Alignment> iter = interval.getAlignmentIterator();
        while(iter.hasNext()){
            Alignment al = iter.next();
            assertNotNull(al);
        }
    }

    /**
     * Load alignment interval. Here for other tests, so we don't need to expose
     * {@link AlignmentDataManager#loadInterval(String, int, int, AlignmentTrack.RenderOptions)}
     */
    public static AlignmentInterval loadInterval(AlignmentDataManager manager, String chr, int start, int end) {
        return manager.loadInterval(chr, start, end, (new AlignmentTrack.RenderOptions()));
    }


    /**
     * Sorts by: read name, start, end
     */
    private class StartEndSorter implements Comparator<Alignment> {

        public int compare(Alignment al1, Alignment al2) {

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
}
