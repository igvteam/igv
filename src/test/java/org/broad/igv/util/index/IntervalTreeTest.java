package org.broad.igv.util.index;

import junit.framework.Assert;
import org.junit.BeforeClass;
import org.junit.Test;


import java.util.List;

/**
 * User: jrobinso
 * Date: Mar 24, 2010
 */
public class IntervalTreeTest {

    static IntervalTree tree;

    @BeforeClass
    public static void setupTree() {
        tree = new IntervalTree();
        tree.insert(new Interval(0, 3, 1));
        tree.insert(new Interval(5, 8, 2));
        tree.insert(new Interval(6, 10, 3));
        tree.insert(new Interval(8, 9, 4));
        tree.insert(new Interval(15, 23, 5));
        tree.insert(new Interval(16, 21, 6));
        tree.insert(new Interval(17, 19, 7));
        tree.insert(new Interval(19, 20, 8));
        tree.insert(new Interval(25, 30, 9));
        tree.insert(new Interval(26, 27, 10));
    }

    @Test
    public void testSearch() {

        List<Interval> intervals = tree.findOverlapping(1, 2);
        Assert.assertNotNull(intervals);
        Assert.assertTrue(intervals.size() > 0);
        for (Interval iv : intervals) {
            Assert.assertTrue(iv.overlaps(1, 2));
        }
    }




}

