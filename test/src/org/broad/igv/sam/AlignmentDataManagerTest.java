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

package org.broad.igv.sam;

import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.Locus;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.sam.reader.ReadGroupFilter;
import org.broad.igv.track.RenderContextImpl;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.io.IOException;
import java.util.*;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012-Jul-12
 */
public class AlignmentDataManagerTest extends AbstractHeadlessTest {

    @Test
    public void testPreloadPanning() throws Exception {

        final String chr = "chr1";
        final int start = 151666494;
        final int halfwidth = 1000;
        final int end = start + 2 * halfwidth;
        int panInterval = halfwidth;

        int numPans = 20 * (end - start) / (panInterval) * AlignmentDataManager.MAX_INTERVAL_MULTIPLE;
        Collection<AlignmentInterval> intervals = AlignmentDataManagerTest.performPanning(chr, start, end, panInterval, numPans);

        assertEquals(1, intervals.size());

    }

    @Test
    public void testPreloadNoMerge() throws Exception {

        final String chr = "chr1";
        final int start = 151666494;
        final int halfwidth = 1000;
        final int end = start + 2 * halfwidth;

        //Load separate intervals, check they don't merge
        AlignmentDataManager manager = getManager171();
        ReferenceFrame frame = new ReferenceFrame("test");
        AlignmentTrack.RenderOptions renderOptions = new AlignmentTrack.RenderOptions();
        frame.setBounds(0, end - start);

        RenderContextImpl context = new RenderContextImpl(null, null, frame, null);

        int lastStart = genome.getChromosome(chr).getLength() - 4 * halfwidth;
        int[] starts = new int[]{500, 5000, 15000, start, 500000, lastStart};
        int[] ends = new int[]{600, 10000, 20000, end, 600000, lastStart + 2 * halfwidth};
        for (int ii = 0; ii < starts.length; ii++) {
            frame.jumpTo(new Locus(chr, starts[ii], ends[ii]));
            int actEnd = (int) frame.getEnd();

            manager.preload(context, renderOptions, false);

            assertManagerHasInterval(manager, chr, starts[ii], actEnd);
        }


    }

    private static void assertManagerHasInterval(AlignmentDataManager manager, String chr, int start, int end) {
        Collection<AlignmentInterval> intervals = manager.getLoadedIntervals();

        boolean haveInterval = false;
        for (AlignmentInterval interval : intervals) {
            haveInterval |= interval.contains(chr, start, end);
        }
        assertTrue(haveInterval);
    }

    public static AlignmentDataManager getManager171() throws IOException {

        String infilepath = TestUtils.LARGE_DATA_DIR + "HG00171.hg18.bam";
        ResourceLocator locator = new ResourceLocator(infilepath);
        AlignmentDataManager manager = new AlignmentDataManager(locator, genome);
        return manager;
    }

    /**
     * Emulates panning across a specific interval.
     *
     * @param chr
     * @param start
     * @param end
     * @param panInterval
     * @param numPans
     * @return
     * @throws IOException
     */
    public static Collection<AlignmentInterval> performPanning(String chr, int start, int end, int panInterval, int numPans) throws IOException {

        AlignmentDataManager manager = getManager171();

        int shift = 0;

        ReferenceFrame frame = new ReferenceFrame("test");
        AlignmentTrack.RenderOptions renderOptions = new AlignmentTrack.RenderOptions();
        frame.setBounds(0, end - start);
        RenderContextImpl context = new RenderContextImpl(null, null, frame, null);

        for (int pp = 0; pp < numPans; pp++) {
            shift = pp * panInterval;
            Locus locus = new Locus(chr, start + shift, end + shift);
            frame.jumpTo(locus);

            manager.preload(context, renderOptions, false);

            assertManagerHasInterval(manager, chr, locus.getStart(), locus.getEnd());
        }

        return manager.getLoadedIntervals();

    }


    @Test
    public void testQuery() throws IOException {
        String testFile = "http://www.broadinstitute.org/igvdata/BodyMap/hg18/50bp/FCA/s_1_1_sequence.bam";
        String sequence = "chr1";
        int start = 44680145;
        int end = 44789983;
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

//            if (i % 2 == 0 && rec.isPaired()) {
//                //Check that paired reads are together
//
//                System.out.println(rec.getReadName());
//
//                System.out.println(result.get(i + 1).getReadName());
//                //assertEquals(rec.getReadName(), result.get(i+1).getReadName());
//            }

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
        PreferenceManager.getInstance().put(PreferenceManager.SAM_MAX_VISIBLE_RANGE, "5");
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

    @Test
    public void testQueryPiledUp() throws Exception {
        PreferenceManager.getInstance().put(PreferenceManager.SAM_MAX_VISIBLE_RANGE, "5");
        PreferenceManager.getInstance().put(PreferenceManager.SAM_DOWNSAMPLE_READS, "false");
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
    @Test
    public void testQueryLargeFile2() throws Exception {
        String path = "http://www.broadinstitute.org/igvdata/1KG/pilot2Bams/NA12878.454.bam";

        ResourceLocator loc = new ResourceLocator(path);

        String sequence = "MT";
        int start = 1000;
        int end = 3000;

        AlignmentDataManager manager = new AlignmentDataManager(loc, genome);

        AlignmentInterval interval = loadInterval(manager, sequence, start, end);

        Iterator<Alignment> iter = interval.getAlignmentIterator();

        int count = 0;
        while (iter.hasNext()) {
            count++;
            Alignment al = iter.next();
            assertTrue(al.getStart() <= end);
            assertTrue(al.getEnd() >= start);
        }

        Assert.assertTrue(count > 0);
        //System.out.println(count + " alignments loaded");

    }

    /**
     * Load alignment interval. Here for other tests, so we don't need to expose
     * {@link AlignmentDataManager#loadInterval(String, int, int, AlignmentTrack.RenderOptions)}
     */
    public static AlignmentInterval loadInterval(AlignmentDataManager manager, String chr, int start, int end) {
        return manager.loadInterval(chr, start, end, new AlignmentTrack.RenderOptions());
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
