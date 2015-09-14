/*
 * The MIT License (MIT)
 *  Copyright (c) 2007-2015 by Institute for Computational Biomedicine,
 *                                          Weill Medical College of Cornell University.
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
