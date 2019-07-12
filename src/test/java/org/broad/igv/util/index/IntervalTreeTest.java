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

