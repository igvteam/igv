/*
 * Copyright (c) 2007-2011 by Institute for Computational Biomedicine,
 *                                          Weill Medical College of Cornell University.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.goby;

import org.broad.igv.feature.LocusScore;
import org.junit.Test;

import java.io.File;
import java.util.List;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * Test that the count goby data source can load data.
 * @author Fabien Campagne
 *         Date: 6/12/11
 *         Time: 1:34 PM
 */
public class GobyCountArchiveDataSourceTest {
    @Test
    public void iteratePickrell() {

        GobyCountArchiveDataSource counts = new GobyCountArchiveDataSource(new File("test/data/goby/counts/GDFQPGI-pickrellNA18486_yale"));
        assertTrue(counts.getAvailableWindowFunctions().size() > 0);
        List<LocusScore> list1 = counts.getSummaryScoresForRange("chr1", 10000, 2000000, 1);
        assertTrue(list1.size() > 0);
        List<LocusScore> list2 = counts.getSummaryScoresForRange("1", 10000, 2000000, 1);
        assertTrue(list2.size() > 0);
        assertEquals(list1.size(), list2.size());
    }

    @Test
    public void iteratePickrellChr3() {

        GobyCountArchiveDataSource counts = new GobyCountArchiveDataSource(new File("test/data/goby/counts/GDFQPGI-pickrellNA18486_yale"));
        assertTrue(counts.getAvailableWindowFunctions().size() > 0);
        List<LocusScore> list1 = counts.getSummaryScoresForRange("chr3", 0,198022431, 0);
        assertTrue(list1.size() > 0);
        List<LocusScore> list2 = counts.getSummaryScoresForRange("3",  0,198022431, 0);
        assertTrue(list2.size() > 0);
        assertEquals(list1.size(), list2.size());
    }
}
