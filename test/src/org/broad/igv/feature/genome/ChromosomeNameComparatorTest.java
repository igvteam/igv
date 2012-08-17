/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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
        int[] expVals = {0, -1, -9, 0, -1, +1, -3, -1, +1};

        ChromosomeNameComparator comp = ChromosomeNameComparator.get();
        for (int ii = 0; ii < set0.length; ii++) {
            assertEquals(expVals[ii], comp.compare(set0[ii], set1[ii]));
        }

        for (int ii = 0; ii < set0.length; ii++) {
            assertEquals(-expVals[ii], comp.compare(set1[ii], set0[ii]));
        }

    }

    @Test
    public void testSort() {

        String[] chrs = {"chr12", "chr10", "chrMT", "chr1", "chrLongName", "chrLongName1"};
        String[] expectedResult = {"chr1", "chr10", "chr12", "chrLongName", "chrLongName1", "chrMT"};

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
