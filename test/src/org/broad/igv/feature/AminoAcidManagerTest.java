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

import org.broad.igv.AbstractHeadlessTest;
import org.junit.Test;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static org.junit.Assert.assertEquals;

/**
 * @author jrobinso
 */
public class AminoAcidManagerTest extends AbstractHeadlessTest {


    public AminoAcidManagerTest() {
    }

    /**
     * Test getting the amino acids for the first exon on gene EGFR
     */
    @Test
    public void testStartEgfr() {
        String expectedSeq = "MRPSGTAGAALLALLAALCPASRALEEKKVC";
        IGVFeature egfr = (IGVFeature) FeatureDB.getFeature("EGFR");

        AminoAcidSequence aaSeq = egfr.getExons().get(0).getAminoAcidSequence(genome);

        int i = 0;
        for (AminoAcid acid : aaSeq.getSequence()) {
            char exp = expectedSeq.charAt(i);
            char actual = acid.getSymbol();
            assertEquals("i=" + i, exp, actual);
            i++;
        }

    }

    @Test
    public void testExon2EGFR() {
        String expectedSeq = "CQGT";
        int expectedOffset = 2;
        IGVFeature egfr = (IGVFeature) FeatureDB.getFeature("EGFR");

        AminoAcidSequence aaSeq = egfr.getExons().get(1).getAminoAcidSequence(genome);
        List<AminoAcid> aaList = aaSeq.getSequence();
        for (int i = 0; i < expectedSeq.length(); i++) {
            assertEquals("i=" + i, expectedSeq.charAt(i), aaList.get(i).getSymbol());

        }

    }

    @Test
    public void testExon3EGFR() {
        final int exonIndex = 2;
        String expectedSeq = "TIQEV";
        int expectedOffset = 0;
        IGVFeature egfr = (IGVFeature) FeatureDB.getFeature("EGFR");

        AminoAcidSequence aaSeq = egfr.getExons().get(exonIndex).getAminoAcidSequence(genome);
        List<AminoAcid> aaList = aaSeq.getSequence();
        for (int i = 0; i < expectedSeq.length(); i++) {
            assertEquals("i=" + i, expectedSeq.charAt(i), aaList.get(i).getSymbol());

        }

    }

    /**
     * Test the last exon of a gene on the negative strand
     */
    @Test
    public void testNegativeEndExon() {
        String expectedSeq = "LLEQNM";
        IGVFeature fbxw7 = (IGVFeature) FeatureDB.getFeature("FBXW7");

        // Insure that this is a negative strand gene, and that its transcript can be read
        assertEquals(Strand.NEGATIVE, fbxw7.getStrand());

        int lastExon = 10;
        AminoAcidSequence aaSeq = fbxw7.getExons().get(lastExon).getAminoAcidSequence(genome);

        int n = expectedSeq.length();
        for (int i = n; i < expectedSeq.length(); i++) {
            assertEquals("# = " + i, expectedSeq.charAt(i - n), aaSeq.getSequence().get(i).getSymbol());
        }
    }

    /**
     * Test a middle exon of a gene on the negative strand (second exon for
     * gene FBXW7)
     */
    @Test
    public void testNegativeGeneMidExon1() {
        String expectedSeq = "QLTQLCQGTKI";
        String expectedEndSeq = "HIGDF";
        IGVFeature fbxw7 = (IGVFeature) FeatureDB.getFeature("FBXW7");


        AminoAcidSequence aaSeq = fbxw7.getExons().get(1).getAminoAcidSequence(genome);


        List<AminoAcid> tmp = aaSeq.getSequence();
        for (int i = 0; i < expectedSeq.length(); i++) {
            assertEquals("# = " + i, expectedSeq.charAt(i), tmp.get(i).getSymbol());
        }

        int delta = tmp.size() - expectedEndSeq.length();
        for (int i = 0; i < expectedEndSeq.length(); i++) {
            assertEquals("# = " + i, expectedEndSeq.charAt(i), tmp.get(i + delta).getSymbol());
        }
    }

    /**
     * Test a middle exon of a gene on the negative strand (third exon for
     * gene FBXW7)
     */
    @Test
    public void testNegativeGeneMidExon2() {
        String expectedSeq = "QLSYV";
        String expectedEndSeq = "GSVV";
        IGVFeature fbxw7 = (IGVFeature) FeatureDB.getFeature("FBXW7");

        AminoAcidSequence aaSeq = fbxw7.getExons().get(2).getAminoAcidSequence(genome);
        List<AminoAcid> tmp = aaSeq.getSequence();
        for (int i = 0; i < expectedSeq.length(); i++) {
            assertEquals("# = " + i, expectedSeq.charAt(i), tmp.get(i).getSymbol());
        }

        int delta = tmp.size() - expectedEndSeq.length();
        for (int i = 0; i < expectedEndSeq.length(); i++) {
            assertEquals("# = " + i, expectedEndSeq.charAt(i), tmp.get(i + delta).getSymbol());
        }

    }

    @Test
    public void translateSequencePositive() {
        // from EGFR start
        String seq = "ATGCGACCC";
        char[] aminoSeq = {'M', 'R', 'P'};

        List<AminoAcid> acids = AminoAcidManager.getAminoAcids(seq, Strand.POSITIVE);
        assertEquals(3, acids.size());
        for (int i = 0; i < 3; i++) {
            assertEquals(aminoSeq[i], acids.get(i).getSymbol());
        }
    }

    @Test
    public void translateSequenceNegative() {
        // From URG$ "end" actually start of coding
        String seq = "CGACGCCAT";
        char[] aminoSeq = {'S', 'A', 'M'};

        List<AminoAcid> acids = AminoAcidManager.getAminoAcids(seq, Strand.NEGATIVE);
        assertEquals(3, acids.size());
        for (int i = 0; i < 3; i++) {
            assertEquals(aminoSeq[i], acids.get(i).getSymbol());
        }
    }

    @Test
    public void testCheckSNPs() {
        Set<String> snps = AminoAcidManager.getAllSNPs("AAA");
        for (String snp : snps) {
            assertEquals(3, snp.length());
            int aas = 0;
            for (char c : snp.toCharArray()) {
                if (c == 'A') {
                    aas += 1;
                }
            }
            assertEquals(2, aas);
        }
    }

    @Test
    public void testgetAminoAcidByName() throws Exception {
        Map<String, String> expected = new HashMap<String, String>();
        expected.put("His", "Histidine");
        expected.put("Thr", "Threonine");
        expected.put("Asn", "Asparagine");
        expected.put("R", "Arginine");
        expected.put("K", "Lysine");
        expected.put("V", "Valine");
        for (String name : expected.keySet()) {
            String exp = expected.get(name);
            String act = AminoAcidManager.getAminoAcidByName(name).getFullName();
            assertEquals(exp, act);
        }
    }
}
