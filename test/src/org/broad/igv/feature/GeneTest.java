/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.feature;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.Feature;
import org.junit.AfterClass;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.BeforeClass;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

/**
 * + example
 * EGFR	NM_005228	chr7	+	55054218	55242524	55054464	55240804	28	55054218,55177472,55178491,55181792,55186480,55187732,55189197,55191016,55191719,55191945,55192849,55195325,55196685,55198919,55200466,55206361,55208169,55209107,55209908,55216479,55226905,55227952,55233903,55235502,55236374,55236921,55237703,55240442,	55054552,55177624,55178675,55181927,55186549,55187851,55189339,55191133,55191846,55192019,55192940,55195525,55196818,55199010,55200624,55206400,55208311,55209230,55210007,55216665,55227061,55228028,55234050,55235600,55236542,55236969,55237812,55242524,
 *
 * @author jrobinso
 */
public class GeneTest {


    static BasicFeature egfr;

    static BasicFeature GTPBP6;
    private static Genome genome;

    public GeneTest() {

    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        genome = TestUtils.loadGenome("hg18");
        egfr = (BasicFeature) FeatureDB.getFeature("egfr");
        GTPBP6 = (BasicFeature) FeatureDB.getFeature("GTPBP6");

    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }


    /**
     * 55054218	55242524	55054464	55240804
     * Test of getStart method, of class IGVFeature.
     */
    @Test
    public void testGetStart() {
        assertEquals(55054218, egfr.getStart());
    }

    /**
     * Test of getCdStart method, of class IGVFeature.
     */
    @Test
    public void testGetCdStart() {
        assertEquals(55054464, egfr.getExons().get(0).getCdStart());
    }

    /**
     * Test of getCdEnd method, of class IGVFeature.
     */
    @Test
    public void testGetCdEnd() {
        assertEquals(55240804, egfr.getExons().get(egfr.getExonCount() - 1).getCdEnd());
    }


    /**
     * Test of getExonCount method, of class IGVFeature.
     */
    @Test
    public void testGetExonCount() {
        assertEquals(28, egfr.getExons().size());
    }

    /**
     * Test of getCodingLength method, of class IGVFeature.
     * <p/>
     * EGFR	NM_005228	chr7	+
     * 55054218	55242525
     * 55054464	55240804	28
     * 55054218,55177472,55178491,55181792,55186480,55187732,55189197,55191016,55191719,55191945,55192849,55195325,55196685,55198919,55200466,55206361,55208169,55209107,55209908,55216479,55226905,55227952,55233903,55235502,55236374,55236921,55237703,55240442,
     * 55054552,55177624,55178675,55181927,55186549,55187851,55189339,55191133,55191846,55192019,55192940,55195525,55196818,55199010,55200624,55206400,55208311,55209230,55210007,55216665,55227061,55228028,55234050,55235600,55236542,55236969,55237812,55242525,
     */
    @Test
    public void testGetCodingLength() {
        int cdStart = 55054464;
        int cdEnd = 55240804;

        int[] exonStarts = {55054218, 55177472, 55178491, 55181792, 55186480,
                55187732, 55189197, 55191016, 55191719, 55191945, 55192849, 55195325,
                55196685, 55198919, 55200466, 55206361, 55208169, 55209107, 55209908,
                55216479, 55226905, 55227952, 55233903, 55235502, 55236374, 55236921,
                55237703, 55240442
        };

        int[] exonEnds = {55054552, 55177624, 55178675, 55181927, 55186549, 55187851,
                55189339, 55191133, 55191846, 55192019, 55192940, 55195525, 55196818, 55199010,
                55200624, 55206400, 55208311, 55209230, 55210007, 55216665, 55227061, 55228028,
                55234050, 55235600, 55236542, 55236969, 55237812, 55242525
        };

        Exon[] exons = new Exon[exonStarts.length];
        for (int i = 0; i < exonStarts.length; i++) {
            Exon exon = new Exon("chr7", exonStarts[i], exonEnds[i], Strand.POSITIVE);
            exon.setCodingStart(cdStart);
            exon.setCodingEnd(cdEnd);
            exons[i] = exon;
        }

        int cdLength0 = exonEnds[0] - cdStart;
        assertEquals(cdLength0, egfr.getExons().get(0).getCodingLength());

        int cdLengthN = cdEnd - exonStarts[exonEnds.length - 1];
        assertEquals(cdLengthN, egfr.getExons().get(egfr.getExonCount() - 1).getCodingLength());

        for (int i = 1; i < egfr.getExonCount() - 1; i++) {
            int expectedLength = exonEnds[i] - exonStarts[i];
            assertEquals(expectedLength, egfr.getExons().get(i).getCodingLength());
        }
    }

    @Test
    public void testGetCodon() {

        NamedFeature feature = FeatureDB.getFeature("egfr");

        // Like nearly all uses of instanceof and casts this is a terrible hack, but one we will live with temporarily
        // until the feature hierarchy can be refactored.
        if (feature instanceof BasicFeature) {

            BasicFeature bf = (BasicFeature) feature;
            int featureStart = bf.getStart();

            // EGFR starts with protiens MRPSG
            Codon codon = bf.getCodon(genome, 2);
            assertEquals('R', codon.getAminoAcid().getSymbol());
            assertEquals(featureStart + 2, codon.getGenomePosition1());
            assertEquals(featureStart + 3, codon.getGenomePosition2());
            assertEquals(featureStart + 4, codon.getGenomePosition3());

            codon = bf.getCodon(genome, 3);
            assertEquals('P', codon.getAminoAcid().getSymbol());
            assertEquals(featureStart + 5, codon.getGenomePosition1());
            assertEquals(featureStart + 6, codon.getGenomePosition2());
            assertEquals(featureStart + 7, codon.getGenomePosition3());

        }

    }


}
