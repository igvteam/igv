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

package org.broad.igv.util.collections;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.data.Interval;
import org.broad.igv.sam.AlignmentDataManager;
import org.broad.igv.sam.AlignmentDataManagerTest;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.junit.After;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012-Oct-16
 */
public class CachedIntervalsTest extends AbstractHeadlessTest {

    private CachedIntervals cachedIntervals;

    @After
    public void tearDown() throws Exception {
        super.tearDown();
        cachedIntervals = null;
    }

    private static List<ReferenceFrame> getTestFrames(String chr) {
        ReferenceFrame frame0 = new ReferenceFrame("leftFrame");
        frame0.setChromosomeName(chr);
        frame0.setOrigin(151666680);
        frame0.setBounds(0, 5);
        System.out.println(frame0.getEnd());

        ReferenceFrame frame1 = new ReferenceFrame("rightFrame");
        frame1.setChromosomeName(chr);
        frame1.setOrigin(153537238);
        frame1.setBounds(0, 5);

        return Arrays.asList(frame0, frame1);
    }

    /**
     * Test that having an interval which spans between two reference frames
     * still exists after trimming
     *
     * @throws Exception
     */
    @Test
    public void testTrimSplit() throws Exception {
        tstTrim(TrimTestType.SPLIT);
    }

    /**
     * Test that having an interval which spans between two reference frames
     * still exists after trimming
     *
     * @throws Exception
     */
    @Test
    public void testTrimSpan() throws Exception {
        tstTrim(TrimTestType.SPAN);
    }

    public void tstTrim(TrimTestType testType) throws Exception {
        //We set a small maxintervalsize to make sure it gets trimmed
        cachedIntervals = new CachedIntervals<Interval>(5, 100);
        String chr = "chr1";

        List<ReferenceFrame> frameList = getTestFrames(chr);
        assertTrue(frameList.get(1).getOrigin() > frameList.get(0).getEnd());
        cachedIntervals.setLocusList(frameList);

        AlignmentDataManager manager = AlignmentDataManagerTest.getManager171();
        int offset = 5;

        int alQueryStart, alQueryEnd;
        int postTrimStart, postTrimEnd;
        switch (testType) {
            case SPLIT:
                alQueryStart = (int) frameList.get(0).getEnd() - 200;
                alQueryEnd = (int) frameList.get(1).getOrigin() + 200;
                postTrimStart = alQueryStart;
                postTrimEnd = alQueryEnd;
                break;
            case SPAN:
                alQueryStart = (int) frameList.get(0).getOrigin() - 100;
                alQueryEnd = (int) frameList.get(1).getEnd() + 100;
                postTrimStart = (int) frameList.get(0).getOrigin();
                postTrimEnd = (int) frameList.get(1).getEnd();
                break;
            default:
                throw new IllegalArgumentException("Unknown testType " + testType);
        }

        //Put in a small interval first, adding the second one triggers the trim
        Interval interval_0 = AlignmentDataManagerTest.loadInterval(manager, chr, alQueryStart, alQueryStart + offset);
        cachedIntervals.put(interval_0);

        Interval interval_1 = AlignmentDataManagerTest.loadInterval(manager, chr, alQueryStart + offset, alQueryEnd + offset);
        cachedIntervals.put(interval_1);

        assertEquals(1, cachedIntervals.getContains(chr, postTrimStart, postTrimEnd, -1).size());
    }

    private enum TrimTestType {
        //Interval is partially covered by two different reference frames
        SPLIT,
        //Interval completely covers two different reference frames
        SPAN
    }
}
