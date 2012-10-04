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

import org.broad.igv.AbstractHeadlessTest;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.*;

/**
 * User: jacob
 * Date: 2012-Oct-04
 */
public class CodonTest extends AbstractHeadlessTest {

    private BasicFeature egfr;

    @Override
    public void setUp() throws Exception {
        super.setUp();
        egfr = (BasicFeature) FeatureDB.getFeature("egfr");
    }

    @Test
    public void testGetCodon() {
        // See http://www.ncbi.nlm.nih.gov/nuccore/NM_201283
        //Note: This covers a break in exons
        char[] expected = "MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYVQRNYDLSFLKTIQEVAGYVLIALNTVERIPLE".toCharArray();
        tstGetCodon(egfr, expected);

    }

    @Test
    public void testGetCodonNeg() {
        //See http://www.ncbi.nlm.nih.gov/nuccore/NM_004985
        String exp_string = "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQV";
        char[] expected = exp_string.toCharArray();

        BasicFeature KRAS = (BasicFeature) FeatureDB.getFeature("KRAS");
        tstGetCodon(KRAS, expected);
    }

    public void tstGetCodon(NamedFeature feature, char[] expected) {
        BasicFeature bf = (BasicFeature) feature;

        for (int pos = 0; pos < expected.length; pos++) {
            Codon codon = bf.getCodon(genome, pos + 1);
            assertEquals(expected[pos], codon.getAminoAcid().getSymbol());
        }
    }

    //Check that there is no exception when loading an
    @Test
    public void testGetCodonInvalid() {

        String geneName = "APOL1";
        int proteinPos = 384;
        //The position isn't valid for this transcript only
        String invalidId = "NM_001136541";
        //valid for these ids:
        //NM_003661, NM_145343, NM_001136540

        List<NamedFeature> featuresList = FeatureDB.getFeaturesList(geneName, 50, false);
        for (NamedFeature feat : featuresList) {
            BasicFeature bf = (BasicFeature) feat;
            System.out.println(bf.getIdentifier());
            Codon c = bf.getCodon(genome, proteinPos);

            if (bf.getIdentifier().equals(invalidId)) {
                assertNull(c);
            } else {
                assertNotNull(c);
                assertTrue(c.isGenomePositionsSet());
            }


        }

    }
}
