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

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.Locus;
import org.broad.igv.track.RenderContextImpl;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Ignore;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static junit.framework.Assert.*;

/**
 * User: jacob
 * Date: 2012-Jul-09
 */
public class AlignmentIntervalTest extends AbstractHeadlessTest {

    @Test
    public void testMergeDense() throws Exception {
        String infilepath = TestUtils.LARGE_DATA_DIR + "HG00171.hg18.bam";
        String chr = "chr1";
        int start = 151666494;
        int halfwidth = 1024;
        tstMerge(infilepath, chr, start, halfwidth);
    }

    @Ignore
    @Test
    public void testMergeSparse() throws Exception {
        String infilepath = "http://www.broadinstitute.org/igvdata/1KG/freeze5_merged/low_coverage_JPT.1.bam";
        String chr = "chr1";
        int halfwidth = 2 * 10000000;
        int start = 151666494 - halfwidth;
        tstMerge(infilepath, chr, start, halfwidth);
    }

    public void tstMerge(String infilepath, String chr, int start, int halfwidth) throws Exception {

        ResourceLocator locator = new ResourceLocator(infilepath);
        AlignmentDataManager baseManager = new AlignmentDataManager(locator, genome);
        AlignmentDataManager testManager = new AlignmentDataManager(locator, genome);

        int end = start + 2 * halfwidth;
        Locus locus = new Locus(chr, start, end);


        ReferenceFrame frame = new ReferenceFrame("test");
        frame.jumpTo(locus);

        frame.setBounds(0, halfwidth);
        RenderContextImpl context = new RenderContextImpl(null, null, frame, null);


        //End is not exact due to zooming
        end = (int) frame.getEnd();

        AlignmentTrack.RenderOptions renderOptions = new AlignmentTrack.RenderOptions();
        baseManager.preload(context, renderOptions, false);

        ArrayList<AlignmentInterval> baseIntervals = (ArrayList) baseManager.getLoadedIntervals();
        assertEquals(1, baseIntervals.size());
        AlignmentInterval baseInterval = baseIntervals.get(0);


        Locus begLocus = new Locus(chr, start, start + halfwidth);
        ReferenceFrame begFrame = new ReferenceFrame(frame);
        begFrame.jumpTo(begLocus);
        RenderContextImpl begContext = new RenderContextImpl(null, null, begFrame, null);

        Locus endLocus = new Locus(chr, start + halfwidth / 2, end);
        ReferenceFrame endFrame = new ReferenceFrame(frame);
        endFrame.jumpTo(endLocus);
        RenderContextImpl endContext = new RenderContextImpl(null, null, endFrame, null);

        testManager.preload(begContext, renderOptions, false);
        ArrayList<AlignmentInterval> begInterval = (ArrayList) testManager.getLoadedIntervals();
        assertEquals(1, begInterval.size());

        testManager.clear();
        testManager.preload(endContext, renderOptions, false);
        ArrayList<AlignmentInterval> endInterval = (ArrayList) testManager.getLoadedIntervals();
        assertEquals(1, endInterval.size());
        AlignmentInterval merged = begInterval.get(0);

        assertTrue(merged.canMerge(endInterval.get(0)));

        assertTrue(merged.merge(endInterval.get(0)));

        assertAlignmentCountsEqual(baseInterval.getCounts(), merged.getCounts());

        TestUtils.assertFeatureListsEqual(baseInterval.getDownsampledIntervals().iterator(), merged.getDownsampledIntervals().iterator());

        TestUtils.assertFeatureListsEqual(baseInterval.getAlignmentIterator(), merged.getAlignmentIterator());

    }

    private void assertAlignmentCountsEqual(AlignmentCounts expected, AlignmentCounts actual) {
        assertEquals(expected.getMaxCount(), actual.getMaxCount());
        assertEquals(expected.getNumberOfPoints(), actual.getNumberOfPoints());
        TestUtils.assertFeaturesEqual(expected, actual);
        byte[] bases = {'a', 'A', 'c', 'C', 'g', 'G', 't', 'T', 'n', 'N'};

        for (int pos = expected.getStart(); pos < expected.getEnd(); pos++) {
            assertEquals(expected.getDelCount(pos), actual.getDelCount(pos));
            assertEquals(expected.getInsCount(pos), actual.getInsCount(pos));

            assertEquals(expected.getTotalCount(pos), actual.getTotalCount(pos));
            assertEquals(expected.getTotalQuality(pos), actual.getTotalQuality(pos));

            //We don't test these because they aren't used anywhere else
            //assertEquals(expected.getPosTotal(pos), actual.getPosTotal(pos));
            //assertEquals(expected.getNegTotal(pos), actual.getNegTotal(pos));

            for (byte b : bases) {
                assertEquals(expected.getCount(pos, b), actual.getCount(pos, b));
                assertEquals(expected.getNegCount(pos, b), actual.getNegCount(pos, b));
                assertEquals(expected.getPosCount(pos, b), actual.getPosCount(pos, b));
                assertEquals(expected.getQuality(pos, b), actual.getQuality(pos, b));

                //Not used anywhere
                //assertEquals(expected.getAvgQuality(pos, b), actual.getAvgQuality(pos, b));
            }
        }

    }

    @Test
    public void testTrimTo() throws Exception {

        String chr = "chr1";
        int start = 151666494;
        int halfwidth = 1000;
        int end = start + 2 * halfwidth;
        int panInterval = halfwidth;
        int buffer = halfwidth / 100;

        int numPans = AlignmentDataManager.MAX_INTERVAL_MULTIPLE;
        List<AlignmentInterval> intervals = (ArrayList<AlignmentInterval>)
                AlignmentDataManagerTest.performPanning(chr, start, end, numPans, panInterval);
        AlignmentInterval interval = intervals.get(0);
        assertTrue(interval.contains(chr, start, end + numPans * panInterval - buffer));


        //Now trim to the middle
        int trimShift = (numPans / 2) * panInterval;
        int trimStart = start + trimShift;
        int trimEnd = end + trimShift;
        assertTrue(interval.trimTo(chr, trimStart, trimEnd, -1));
        assertFalse(interval.trimTo(chr, trimStart, trimEnd, -1));
        assertFalse(interval.trimTo(chr, trimStart, trimEnd, -1));
        assertFalse(interval.trimTo(chr, trimStart, trimEnd, -1));

        assertFalse(interval.contains(chr, start, end + numPans * panInterval));
        assertTrue(interval.contains(chr, trimStart, trimEnd, -1));

    }


}
