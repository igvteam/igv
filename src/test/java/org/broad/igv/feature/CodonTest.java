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
