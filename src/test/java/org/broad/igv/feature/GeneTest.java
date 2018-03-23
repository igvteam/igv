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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.feature;

import org.broad.igv.AbstractHeadlessTest;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * + example
 * EGFR	NM_005228	chr7	+	55054218	55242524	55054464	55240804	28	55054218,55177472,55178491,55181792,55186480,55187732,55189197,55191016,55191719,55191945,55192849,55195325,55196685,55198919,55200466,55206361,55208169,55209107,55209908,55216479,55226905,55227952,55233903,55235502,55236374,55236921,55237703,55240442,	55054552,55177624,55178675,55181927,55186549,55187851,55189339,55191133,55191846,55192019,55192940,55195525,55196818,55199010,55200624,55206400,55208311,55209230,55210007,55216665,55227061,55228028,55234050,55235600,55236542,55236969,55237812,55242524,
 *
 * @author jrobinso
 */
public class GeneTest extends AbstractHeadlessTest {


    static BasicFeature egfr;

    static BasicFeature GTPBP6;

    public GeneTest() {

    }

    public static void main(String[] args) {
        //        int[] positions = new int[]{0,1,2,3,4,5,6,7};
//
//        String ret = "";
//        for(int pos: positions){
//            Codon codon = bf.getCodon(genome, pos + 1);
//            //System.out.println("Pos: " + (pos + 1) + " AA: " + codon.getAminoAcid().getSymbol());
//            ret += codon.getAminoAcid().getSymbol();
//        }
//        System.out.println(ret);
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        AbstractHeadlessTest.setUpClass();
        egfr = (BasicFeature) FeatureDB.getFeature("egfr");
        GTPBP6 = (BasicFeature) FeatureDB.getFeature("GTPBP6");

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

//        Exon[] exons = new Exon[exonStarts.length];
//        for (int i = 0; i < exonStarts.length; i++) {
//            Exon exon = new Exon("chr7", exonStarts[i], exonEnds[i], Strand.POSITIVE);
//            exon.setCodingStart(cdStart);
//            exon.setCodingEnd(cdEnd);
//            exons[i] = exon;
//        }

        int cdLength0 = exonEnds[0] - cdStart;
        assertEquals(cdLength0, egfr.getExons().get(0).getCodingLength());

        int cdLengthN = cdEnd - exonStarts[exonEnds.length - 1];
        assertEquals(cdLengthN, egfr.getExons().get(egfr.getExonCount() - 1).getCodingLength());

        for (int i = 1; i < egfr.getExonCount() - 1; i++) {
            int expectedLength = exonEnds[i] - exonStarts[i];
            assertEquals(expectedLength, egfr.getExons().get(i).getCodingLength());
        }
    }


}
