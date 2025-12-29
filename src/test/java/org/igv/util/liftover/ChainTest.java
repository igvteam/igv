package org.igv.util.liftover;

import org.igv.feature.Range;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import static org.junit.Assert.*;

public class ChainTest {


    /**
     * Span contained entirely in block
     */
    @Test
    public void contains() {

        Chain chain = getTestChain();

        Range span = new Range("chr1",4784516, 4784516 + 10);
        List<Range> mapped = chain.map(span);

        assertEquals(1, mapped.size());
        Range m = mapped.get(0);

        assertEquals(0, m.start);
        assertEquals(10, (m.end - m.start));
    }

    /**
     * Span overlaps end of first interval by +/- 10
     */
    @Test
    public void overlappingEnd() {

        Chain chain = getTestChain();

        int boundary = 4784516 + 1725;
        Range span = new Range("chr1", boundary - 10, boundary + 10);
        List<Range> mapped = chain.map(span);

        assertEquals(1, mapped.size());
        Range m = mapped.get(0);

        assertEquals(1725 - 10, m.start);
        assertEquals(10, (m.end - m.start));
    }

    /**
     * Span overlaps start of first interval by +/- 10
     */
    @Test
    public void overlappingStart() {
        Chain chain = getTestChain();

        int boundary = 4784516;
        Range span = new Range("chr1",boundary - 10, boundary + 10);
        List<Range> mapped = chain.map(span);

        assertEquals(1, mapped.size());
        Range m = mapped.get(0);

        assertEquals(0, m.start);
        assertEquals(10, (m.end - m.start));
    }

    /**
     * Overlaps 2 intervals
     */
    @Test
    public void overlappingTwoIntervals() {
        Chain chain = getTestChain();

        int boundary = 4784516 + 1725 + 14679;
        Range span = new Range("chr1",4784516 - 10, boundary + 10);
        List<Range> mapped = chain.map(span);

        assertEquals(2, mapped.size());

        Collections.sort(mapped, (o1, o2) -> o1.start - o2.start);

        Range m1 = mapped.get(0);
        assertEquals(0, m1.start);
        assertEquals(1725, (m1.end - m1.start));

        Range m2 = mapped.get(1);
        assertEquals(1725, m2.start);
        assertEquals(10, (m2.end - m2.start));
    }

    /**
     * Outside all intervals
     */
    @Test
    public void outside() {
        Chain chain = getTestChain();

        Range span = new Range("chr1",10,20);
        List<Range> mapped = chain.map(span);
        assertEquals(0, mapped.size());
    }



    private Chain getTestChain() {

        Chain chain = new Chain("chr1", 190687455, 4784516, 190687455, "chr1", 1650847, 0, 1650847, "1");

        List<String []> alignments = new ArrayList<>();
        alignments.add(new String [] {"1725", "14679", "0"});
        alignments.add(new String [] {"3914", "2468", "0"});
        alignments.add(new String [] {"1535"});

        chain.setAlignments(alignments);

        return chain;

    }
}