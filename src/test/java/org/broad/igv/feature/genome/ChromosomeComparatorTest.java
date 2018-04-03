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

package org.broad.igv.feature.genome;

import org.broad.igv.feature.Chromosome;
import org.junit.Test;

import static junit.framework.Assert.assertTrue;

/**
 * @author jrobinso
 *         Date: 2/13/13
 *         Time: 2:21 PM
 */
public class ChromosomeComparatorTest {


    // Test chr names which resolve to numbers too large for "int".  In response to bug report.

    @Test
    public void testCompareLargeNumbers() throws Exception {
        long num1 = 7000000037415152l;
        long num2 = 7000000037415153l;
        ChromosomeComparator comp = new ChromosomeComparator(0);
        Chromosome chr1 = new Chromosome(0, String.valueOf(num1), 1 );
        Chromosome chr2 = new Chromosome(0, String.valueOf(num2), 1 );
        int value = comp.compare(chr1, chr2);
        assertTrue(value < 0);
    }
}
