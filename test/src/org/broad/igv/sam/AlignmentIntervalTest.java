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
import org.junit.Test;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import static junit.framework.Assert.*;

/**
 * User: jacob
 * Date: 2012-Jul-09
 */
public class AlignmentIntervalTest extends AbstractHeadlessTest {

    @Test
    public void testMerge() throws Exception {
        String infilepath = TestUtils.LARGE_DATA_DIR + "HG00171.hg18.bam";
        ResourceLocator locator = new ResourceLocator(infilepath);
        AlignmentDataManager baseManager = new AlignmentDataManager(locator, genome);
        AlignmentDataManager testManager = new AlignmentDataManager(locator, genome);

        String chr = "chr1";
        int start = 151666494;
        int halfwidth = 1000;
        int end = start + 2 * halfwidth;
        Locus locus = new Locus(chr, start, end);


        ReferenceFrame frame = new ReferenceFrame("test");
        frame.setInterval(locus);

        frame.setBounds(0, halfwidth);
        RenderContextImpl context = new RenderContextImpl(null, null, frame, null);

        AlignmentTrack.RenderOptions renderOptions = new AlignmentTrack.RenderOptions();
        baseManager.preload(context, renderOptions, null, false);

        ArrayList<AlignmentInterval> baseIntervals = (ArrayList) baseManager.getLoadedIntervals();
        assertEquals(1, baseIntervals.size());
        AlignmentInterval base = baseIntervals.get(0);


        Locus begLocus = new Locus(chr, start, start + halfwidth);
        ReferenceFrame begFrame = new ReferenceFrame(frame);
        begFrame.setInterval(begLocus);
        RenderContextImpl begContext = new RenderContextImpl(null, null, begFrame, null);

        Locus endLocus = new Locus(chr, start + halfwidth / 2, end);
        ReferenceFrame endFrame = new ReferenceFrame(frame);
        endFrame.setInterval(endLocus);
        RenderContextImpl endContext = new RenderContextImpl(null, null, endFrame, null);

        testManager.preload(begContext, renderOptions, null, false);
        ArrayList<AlignmentInterval> begInterval = (ArrayList) testManager.getLoadedIntervals();
        assertEquals(1, begInterval.size());

        testManager.clear();
        testManager.preload(endContext, renderOptions, null, false);
        ArrayList<AlignmentInterval> endInterval = (ArrayList) testManager.getLoadedIntervals();
        assertEquals(1, endInterval.size());
        AlignmentInterval merged = begInterval.get(0);
        merged.merge(endInterval.get(0));

        assertEquals(base.getCounts().size(), merged.getCounts().size());

        Iterator<Alignment> iter1 = baseIntervals.get(0).getAlignmentIterator();
        Iterator<Alignment> iter2 = merged.getAlignmentIterator();

        Alignment expected, actual;
        int count = 0;
        //Because we sort by start position, and some alignments have the same start position,
        //there can be some ambiguitiy in which comes first, so we don't compare directly.
        //Just make sure they end up the same in the end.
        Set<String> expectedReadNames = new HashSet<String>();
        Set<String> actualReadNames = new HashSet<String>();
        while (iter1.hasNext()) {
            assertTrue("Not enough alignments, empty at count " + count, iter2.hasNext());
            expected = iter1.next();
            actual = iter2.next();

            assertEquals(expected.getAlignmentStart(), actual.getAlignmentStart());
            assertEquals(expected.getChr(), actual.getChr());
            count++;

            expectedReadNames.add(expected.getReadName());
            actualReadNames.add(actual.getReadName());
        }

        //System.out.println(count + " Alignments");
        assertFalse(iter2.hasNext());
        assertTrue("No data loaded", count > 0);

        assertEquals(expectedReadNames, actualReadNames);
    }
}
