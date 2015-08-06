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
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static org.junit.Assert.*;

/**
 * @author jrobinso
 */
public class AminoAcidManagerTest extends AbstractHeadlessTest {


    @Override
    public void setUp() throws Exception {
        super.setUp();
        AminoAcidManager.resetToDefaultCodonTables();
    }

    /**
     * Test getting the amino acids for the first exon on gene EGFR
     */
    @Test
    public void testStartEgfr() {
        String expectedSeq = "MRPSGTAGAALLALLAALCPASRALEEKKVC";
        IGVFeature egfr = (IGVFeature) FeatureDB.getFeature("EGFR");

        Exon nextExon = egfr.getExons().get(1);
        AminoAcidSequence aaSeq = egfr.getExons().get(0).getAminoAcidSequence(genome, null, nextExon);

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
        // Note:  readingFrame == 2
        String expectedSeq = "VCQGT";
        IGVFeature egfr = (IGVFeature) FeatureDB.getFeature("EGFR");

        Exon exon = egfr.getExons().get(1);
        Exon prevExon = egfr.getExons().get(0);
        Exon nextExon = egfr.getExons().get(2);

        AminoAcidSequence aaSeq = exon.getAminoAcidSequence(genome, prevExon, nextExon);
        List<AminoAcid> aaList = aaSeq.getSequence();
        for (int i = 0; i < expectedSeq.length(); i++) {
            assertEquals("i=" + i, expectedSeq.charAt(i), aaList.get(i).getSymbol());
        }

    }

    @Test
    public void testExon3EGFR() {
        final int exonIndex = 2;
        String expectedSeq = "TIQEV";
        String expectedEndSeq = "NLQE";

        IGVFeature egfr = (IGVFeature) FeatureDB.getFeature("EGFR");
        Exon exon = egfr.getExons().get(exonIndex);
        Exon prevExon = egfr.getExons().get(exonIndex-1);
        Exon nextExon = egfr.getExons().get(exonIndex+1);

        AminoAcidSequence aaSeq = exon.getAminoAcidSequence(genome, prevExon, nextExon);
        List<AminoAcid> aaList = aaSeq.getSequence();
        for (int i = 0; i < expectedSeq.length(); i++) {
            assertEquals("i=" + i, expectedSeq.charAt(i), aaList.get(i).getSymbol());

        }
        for (int i = 0; i < expectedEndSeq.length(); i++) {
            int aaIdx = aaList.size() - expectedEndSeq.length() + i;
            assertEquals("i=" + i, expectedEndSeq.charAt(i), aaList.get(aaIdx).getSymbol());

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

        int lastExonIdx = 10;
        Exon exon = fbxw7.getExons().get(lastExonIdx);
        Exon prevExon = fbxw7.getExons().get(lastExonIdx-1);
        Exon nextExon = null;


        AminoAcidSequence aaSeq = exon.getAminoAcidSequence(genome, prevExon, nextExon);

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

        String expectedSeq = "GQLTQLCQGTKI";
        String expectedEndSeq = "HIGDF";

        IGVFeature fbxw7 = (IGVFeature) FeatureDB.getFeature("FBXW7");
        Exon exon = fbxw7.getExons().get(1);
        Exon prevExon = fbxw7.getExons().get(0);
        Exon nextExon = fbxw7.getExons().get(2);


        AminoAcidSequence aaSeq = exon.getAminoAcidSequence(genome, prevExon, nextExon);


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
        String expectedEndSeq = "GSVVR";
        IGVFeature fbxw7 = (IGVFeature) FeatureDB.getFeature("FBXW7");

        Exon exon = fbxw7.getExons().get(2);
        Exon prevExon = fbxw7.getExons().get(1);
        Exon nextExon = fbxw7.getExons().get(3);

        AminoAcidSequence aaSeq = exon.getAminoAcidSequence(genome, prevExon, nextExon);
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

        List<AminoAcid> acids = AminoAcidManager.getInstance().getAminoAcids(Strand.POSITIVE, seq);
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

        List<AminoAcid> acids = AminoAcidManager.getInstance().getAminoAcids(Strand.NEGATIVE, seq);
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

    //Just check that we have the ids we expect
    @Test
    public void testCodonTablesExist() throws Exception {
        int[] expIds = {1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24};
        for (int id : expIds) {
            assertTrue(AminoAcidManager.getInstance().setCodonTable(AminoAcidManager.DEFAULT_CODON_TABLE_PATH, id));
        }
    }

    //Sample a few codon tables, test different translations
    @Test
    public void testVariousCodonTables() throws Exception {

        int[] codonTableIds = {2, 2, 3, 6, 16, 22, 23, 24};
        String[] testCodons = {"AGA", "TGA", "CTG", "TAA", "TAG", "TCA", "TTA", "AGA"};
        char[] expAAs = {'*', 'W', 'T', 'Q', 'L', '*', '*', 'S'};

        for (int ii = 0; ii < expAAs.length; ii++) {
            AminoAcidManager aam = AminoAcidManager.getInstance();
            int id = codonTableIds[ii];
            boolean loaded = aam.setCodonTable(AminoAcidManager.DEFAULT_CODON_TABLE_PATH, id);

            assertTrue("Failed to load codon table with id " + id, loaded);

            AminoAcid actualAA = aam.getAminoAcid(testCodons[ii]);
            assertNotSame("Got null amino acid for " + testCodons[ii], AminoAcid.NULL_AMINO_ACID, actualAA);
            assertEquals(String.valueOf(expAAs[ii]), String.valueOf(actualAA.getSymbol()));

            //We want to only store one copy of each amino acid. Check that this is the case
            assertTrue(AminoAcidManager.getAminoAcidByName(actualAA.getShortName()) == actualAA);
        }
    }

}
