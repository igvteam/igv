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

import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * User: jacob
 * Date: 2012-Aug-17
 */
public class ChromosomeNameComparatorTest {

    @Test
    public void testCompare() throws Exception {
        String[] set0 = {"chr1", "chr1", "chr1", "CHR1", "chr2a", "chrM", "chrUn_12_gl129"};
        String[] set1 = {"chr1", "chr2", "chr10", "chr1", "chr3", "taoheu", "chrXn_01"};
        int[] expVals = {0, -1, -9, 0, -1, +1, -3};

        ChromosomeNameComparator comp = ChromosomeNameComparator.get();
        for (int ii = 0; ii < set0.length; ii++) {
            assertEquals(expVals[ii], comp.compare(set0[ii], set1[ii]));
        }

        for (int ii = 0; ii < set0.length; ii++) {
            assertEquals(-expVals[ii], comp.compare(set1[ii], set0[ii]));
        }

    }


}
