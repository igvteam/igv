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
    public void testSor1() throws Exception {
        String[] chrs = {"chr1", "chr2a", "chrM", "chrUn_12_gl129", "scaf1_100_b12", "scaf2_200_b80"};
        String[] expectedResult = {"chr1", "chr2a", "chrM", "chrUn_12_gl129", "scaf1_100_b12", "scaf2_200_b80"};
        Arrays.sort(chrs, ChromosomeNameComparator.get());
        for (int i = 0; i < chrs.length; i++) {
            Assert.assertEquals(expectedResult[i], chrs[i]);
        }

    }


    @Test
    public void foo() {
        ChromosomeNameComparator comparator = ChromosomeNameComparator.get();
        System.out.println(comparator.compare("chrUr", "scaffold"));
        System.out.println(comparator.compare("MT", "scaffold"));
    }

    @Test
    public void testSort2() {

        String[] chrs = {"chr12", "chr10", "chr2b", "chr2a", "chrX", "chrM", "chr1", "chrLongName", "chrLongName1", "scaffold"};
        String[] expectedResult = {"chr1", "chr2a", "chr2b", "chr10", "chr12", "chrX", "chrM", "chrLongName", "chrLongName1", "scaffold"};

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

    @Test
    public void testSort3() {

        String[] chrs = {"1", "10", "2b", "2a", "MT", "scaffold"};
        String[] expectedResult = {"1", "2a", "2b", "10", "MT", "scaffold"};

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
