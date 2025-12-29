package org.igv.util.liftover;

import org.igv.feature.Range;
import org.igv.util.TestUtils;
import org.junit.Test;

import java.io.IOException;
import java.util.Collections;
import java.util.List;

import static org.junit.Assert.*;

/**
 * Test Liftover class.  More extensive tests of mapping are in ChainTest, the purpose of this test is to verify the
 * liftover file was parsed correctly.
 *
 */
public class LiftoverTest {

/*
chain 1000 chr1 190687455 + 4784516 190687455 chr1 1650847 + 0 1650847 1
1725 14679 0
3914 2468 0
1535 14340 0
501
*/

    @Test
    public void testLiftover1() throws IOException {

        String liftoverFile = TestUtils.DATA_DIR + "liftover/mm10_K27Ac.chain";
        Liftover liftover = Liftover.load(liftoverFile);

        Range span = new Range("chr1",4784516, 4784516 + 10);
        List<Range> mapped = liftover.map(span);

        assertEquals(1, mapped.size());
        Range m = mapped.get(0);

        assertEquals(0, m.start);
        assertEquals(10, (m.end - m.start));

    }

    @Test
    public void testLiftover2() throws IOException {

/*
chain 1000 chr12 117027732 + 3101290 117027732 chr12 930390 + 0 930390 4
1483 132410 0
1947 1244 0
630 70035 0
1781 5972 0
501
*/

        String liftoverFile = TestUtils.DATA_DIR + "liftover/mm10_K27Ac.chain";
        Liftover liftover = Liftover.load(liftoverFile);

        Range span = new Range("chr12",3101290, 3101290 + 10);
        List<Range> mapped = liftover.map(span);

        assertEquals(1, mapped.size());
        Range m = mapped.get(0);

        assertEquals(0, m.start);
        assertEquals(10, (m.end - m.start));

    }

    /**
     * Test feature that overlaps 2 intervals, spanning the entirety of the first interval and first 10 bases of second.
     * As the mapped regions are contiguous they should be combined to yield one range with size
     *     (total size of first interval) + (10 bases of second)
     * @throws IOException
     */
    @Test
    public void overlappingTwoIntervals() throws IOException {

        String liftoverFile = TestUtils.DATA_DIR + "liftover/mm10_K27Ac.chain";
        Liftover liftover = Liftover.load(liftoverFile);

        final int firstTargetStart = 4784516;
        int secondTargetStart = firstTargetStart + 1725 + 14679;
        Range span = new Range("chr1", firstTargetStart - 10, secondTargetStart + 10);
        List<Range> mapped = liftover.map(span);

        assertEquals(1, mapped.size());

        Range m1 = mapped.get(0);
        int expectedSize = 1725 + 10;
        assertEquals(0, m1.start);
        assertEquals(expectedSize, (m1.end - m1.start));

    }

}