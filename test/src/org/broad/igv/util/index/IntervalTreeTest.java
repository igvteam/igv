/*
 * Copyright (c) 2007-2010 by The Broad Institute, Inc. and the Massachusetts Institute of Technology.
 * All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL), Version 2.1 which
 * is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR WARRANTIES OF
 * ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT
 * OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR
 * RESPECTIVE TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES OF
 * ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES, ECONOMIC
 * DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER THE BROAD OR MIT SHALL
 * BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE
 * FOREGOING.
 */

package org.broad.igv.util.index;

import junit.framework.Assert;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.CloseableTribbleIterator;
import org.broad.tribble.FeatureReader;
import org.broad.tribble.bed.BEDCodec;
import org.broad.tribble.bed.BEDFeature;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.junit.BeforeClass;
import org.junit.Test;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static org.junit.Assert.assertTrue;

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

