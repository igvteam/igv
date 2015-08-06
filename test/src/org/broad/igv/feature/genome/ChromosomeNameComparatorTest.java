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

import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;

import static junit.framework.Assert.assertEquals;

/**
 * User: jacob
 * Date: 2012-Aug-17
 */
public class ChromosomeNameComparatorTest {

    @Test
    public void testCompare() throws Exception {
        String[] set0 = {"chr1", "chr1", "chr1", "CHR1", "chr2a", "chrM", "chrUn_12_gl129", "scaf1_100_b12",
                "scaf2_200_b80"};
        String[] set1 = {"chr1", "chr2", "chr10", "chr1", "chr3", "taoheu", "chrXn_01", "scaf1_101_b50",
                "scaf1_100_b12"};
        int[] expVals = {0, -1, -1, 0, -1, +1, -1, -1, +1};

        ChromosomeNameComparator comp = ChromosomeNameComparator.get();
        for (int ii = 0; ii < set0.length; ii++) {
            // Compare signs
            int tmp = comp.compare(set0[ii], set1[ii]) ;
            int sign = (tmp == 0 ? 0 : (tmp < 0 ? -1 : 1));
            assertEquals(expVals[ii], sign);
        }

        for (int ii = 0; ii < set0.length; ii++) {
            int tmp = comp.compare(set1[ii], set0[ii]) ;
            int sign = (tmp == 0 ? 0 : (tmp < 0 ? -1 : 1));
            assertEquals(-expVals[ii], sign);
        }

    }

    @Test
    public void testSort() {

        String[] chrs = {"chr12", "chr10", "chr2", "chrX", "chrM", "chr1", "chrLongName", "chrLongName1"};
        String[] expectedResult = {"chr1", "chr2", "chr10", "chr12", "chrLongName", "chrLongName1", "chrX", "chrM"};

        Arrays.sort(chrs, ChromosomeNameComparator.get());
        for (int i = 0; i < chrs.length; i++) {
            Assert.assertEquals(expectedResult[i], chrs[i]);
        }
        System.out.println();

        chrs = new String[]{"scaffold_v2_10414", "scaffold_v2_100", "scaffold_v2_101", "scaffold_v2_10415"};
        expectedResult = new String[]{"scaffold_v2_100", "scaffold_v2_101", "scaffold_v2_10414", "scaffold_v2_10415"};
        Arrays.sort(chrs, ChromosomeNameComparator.get());
        for (int i = 0; i < chrs.length; i++) {
            Assert.assertEquals(expectedResult[i], chrs[i]);
        }

    }


}
