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
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.tools.IgvTools;
import org.broad.igv.track.RenderContextImpl;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.IOException;
import java.util.Collection;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012-Jul-12
 */
public class AlignmentDataManagerTest extends AbstractHeadlessTest {

    @Test
    public void testPreload() throws Exception {
        String chr = "chr1";
        int start = 151666494;
        int halfwidth = 1000;
        int end = start + 2 * halfwidth;
        int panInterval = halfwidth;

        int numPans = 20 * (end - start) / (panInterval) * AlignmentDataManager.MAX_INTERVAL_MULTIPLE;
        Collection<AlignmentInterval> intervals = AlignmentDataManagerTest.performPanning(chr, start, end, panInterval, numPans);

        assertEquals(1, intervals.size());

        //Load separate intervals, check they don't merge
        AlignmentDataManager manager = getManager171();
        ReferenceFrame frame = new ReferenceFrame("test");
        AlignmentTrack.RenderOptions renderOptions = new AlignmentTrack.RenderOptions();
        frame.setBounds(0, end - start);
        RenderContextImpl context = new RenderContextImpl(null, null, frame, null);

        int[] starts = new int[]{500, 5000, 15000, start, 500000, start * 2};
        int[] ends = new int[]{600, 10000, 20000, end, 600000, start * 2 + halfwidth};
        for (int ii = 0; ii < starts.length; ii++) {
            frame.setInterval(new Locus(chr, starts[ii], ends[ii]));
            manager.preload(context, renderOptions, null, false);

            assertManagerHasInterval(manager, chr, starts[ii], ends[ii]);
        }


    }

    private static void assertManagerHasInterval(AlignmentDataManager manager, String chr, int start, int end) {
        Collection<AlignmentInterval> intervals = manager.getLoadedIntervals();
        assertTrue(intervals.size() < AlignmentDataManager.CACHE_SIZE);

        boolean haveInterval = false;
        for (AlignmentInterval interval : intervals) {
            haveInterval |= interval.contains(chr, start, end);
        }
        assertTrue(haveInterval);
    }

    private static AlignmentDataManager getManager171() throws IOException {

        String infilepath = TestUtils.LARGE_DATA_DIR + "HG00171.hg18.bam";
        Genome genome = IgvTools.loadGenome("hg18", true);
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
            frame.setInterval(locus);

            manager.preload(context, renderOptions, null, false);

            assertManagerHasInterval(manager, chr, locus.getStart(), locus.getEnd());
        }

        return manager.getLoadedIntervals();

    }
}
