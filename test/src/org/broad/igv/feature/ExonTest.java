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

package org.broad.igv.feature;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.TestUtils;
import org.junit.After;
import org.junit.BeforeClass;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * User: jacob
 * Date: 2012/05/17
 */
public class ExonTest {

    static Genome genome;

    @BeforeClass
    public static void setUp() throws Exception {
        TestUtils.setUpHeadless();
        genome = TestUtils.loadGenome();

    }

    @After
    public void tearDown() throws Exception {

    }

    @Test
    public void testGetAminoAcidNumberPos() throws Exception {
        tstGetAminoAcidNumber("EGFR", 0, true);
    }

    @Test
    public void testGetAminoAcidNumberNeg() throws Exception {
        tstGetAminoAcidNumber("KRAS", -2, false);
    }

    /**
     * exonNum is interpreted as counting backwards
     * if negative. So -1 is the last one.
     *
     * @param gene
     * @param exonNum
     * @param positive
     * @throws Exception
     */
    public void tstGetAminoAcidNumber(String gene, int exonNum, boolean positive) throws Exception {

        int[] genomicOffsets = new int[]{0, 1, 2, 3, 4, 5};
        int[] expAANumbers = new int[]{1, 1, 1, 2, 2, 2};

        BasicFeature egfr = (BasicFeature) FeatureDB.getFeature(gene);
        if (exonNum < 0) {
            exonNum = egfr.getExonCount() + exonNum;
        }
        Exon exon = egfr.getExons().get(exonNum);
        int ind = 0;
        int mult = positive ? 1 : -1;
        int start = positive ? exon.getCdStart() : exon.getCdEnd();
        for (Integer offset : genomicOffsets) {
            int genomicPosition = start + mult * offset;
            int AANumber = exon.getAminoAcidNumber(genomicPosition);
            assertEquals(expAANumbers[ind], AANumber);
            ind++;

        }

    }

}
