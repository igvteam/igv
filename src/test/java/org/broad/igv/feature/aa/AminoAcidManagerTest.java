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
package org.broad.igv.feature.aa;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.*;
import org.broad.igv.feature.tribble.UCSCGeneTableCodec;
import org.junit.Assert;
import org.junit.Ignore;
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
        CodonTableManager.getInstance().resetToDefaults();
    }

    @Test
    public void testFoo() {

        String expectedSeq = "LRAAK";

        UCSCGeneTableCodec codec = new UCSCGeneTableCodec(UCSCGeneTableCodec.Type.GENEPRED, null);
        String str = "1766\tNM_182679\tchr1\t-\t154830723\t154837903\t154831628\t154835462\t8\t154830723,154832496,154832802,154834450,154834642,154834836,154835383,154837824,\t154832265,154832547,154832894,154834537,154834735,154834910,154835520,154837903,\t0\tGPATCH4\tcmpl\tcmpl\t2,2,0,0,0,1,0,-1,";
        BasicFeature feature = codec.decode(str);

        Exon exon = feature.getExons().get(0);
        Exon prevExon = null;
        Exon nextExon = feature.getExons().get(1);
        AminoAcidSequence aaSeq = exon.getAminoAcidSequence(genome, prevExon, nextExon);
        assertNotNull(aaSeq);

        int aaSeqLen = aaSeq.getSequence().size();
        int start = aaSeqLen - expectedSeq.length();
        for (int i = start; i < aaSeqLen - 1; i++) {
            char exp = expectedSeq.charAt(i - start);
            char actual = aaSeq.getSequence().get(i).getAminoAcid().getSymbol();
            assertEquals("i=" + i, exp, actual);
        }
    }

    /**
     * Test getting the amino acids for the first exon on gene EGFR
     */
    @Test
    public void testStartEgfr() {
        String expectedSeq = "MRPSGTAGAALLALLAALCPASRALEEKKVC";
        IGVFeature egfr = (IGVFeature) genome.getFeatureDB().getFeature("EGFR");

        Exon nextExon = egfr.getExons().get(1);
        AminoAcidSequence aaSeq = egfr.getExons().get(0).getAminoAcidSequence(genome, null, nextExon);

        int i = 0;
        for (CodonAA acid : aaSeq.getSequence()) {
            char exp = expectedSeq.charAt(i);
            char actual = acid.getSymbol();
            assertEquals("i=" + i, exp, actual);
            i++;
        }

    }

    @Test
    public void testExon2EGFR() {
        // Note:  frame == 1,  phase == 2
        String expectedSeq = "VCQGT";
        IGVFeature egfr = (IGVFeature) genome.getFeatureDB().getFeature("EGFR");

        Exon exon = egfr.getExons().get(1);
        Exon prevExon = egfr.getExons().get(0);
        Exon nextExon = egfr.getExons().get(2);

        AminoAcidSequence aaSeq = exon.getAminoAcidSequence(genome, prevExon, nextExon);
        List<CodonAA> aaList = aaSeq.getSequence();
        for (int i = 0; i < expectedSeq.length(); i++) {
            assertEquals("i=" + i, expectedSeq.charAt(i), aaList.get(i).getSymbol());
        }

    }

    @Test
    public void testExon3EGFR() {
        final int exonIndex = 2;
        String expectedSeq = "TIQEV";
        String expectedEndSeq = "NLQE";

        IGVFeature egfr = (IGVFeature) genome.getFeatureDB().getFeature("EGFR");
        Exon exon = egfr.getExons().get(exonIndex);
        Exon prevExon = egfr.getExons().get(exonIndex - 1);
        Exon nextExon = egfr.getExons().get(exonIndex + 1);

        AminoAcidSequence aaSeq = exon.getAminoAcidSequence(genome, prevExon, nextExon);
        List<CodonAA> aaList = aaSeq.getSequence();
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
        IGVFeature fbxw7 = (IGVFeature) genome.getFeatureDB().getFeature("FBXW7");

        // Insure that this is a negative strand gene, and that its transcript can be read
        Assert.assertEquals(Strand.NEGATIVE, fbxw7.getStrand());

        int lastExonIdx = 10;
        Exon exon = fbxw7.getExons().get(lastExonIdx);
        Exon prevExon = fbxw7.getExons().get(lastExonIdx - 1);
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

        IGVFeature fbxw7 = (IGVFeature) genome.getFeatureDB().getFeature("FBXW7");
        Exon exon = fbxw7.getExons().get(1);
        Exon prevExon = fbxw7.getExons().get(0);
        Exon nextExon = fbxw7.getExons().get(2);


        AminoAcidSequence aaSeq = exon.getAminoAcidSequence(genome, prevExon, nextExon);


        List<CodonAA> tmp = aaSeq.getSequence();
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
        IGVFeature fbxw7 = (IGVFeature) genome.getFeatureDB().getFeature("FBXW7");

        Exon exon = fbxw7.getExons().get(2);
        Exon prevExon = fbxw7.getExons().get(1);
        Exon nextExon = fbxw7.getExons().get(3);

        AminoAcidSequence aaSeq = exon.getAminoAcidSequence(genome, prevExon, nextExon);
        List<CodonAA> tmp = aaSeq.getSequence();
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

        CodonTable defaultTable = CodonTableManager.getInstance().getDefaultCodonTable();
        List<CodonAA> acids = AminoAcidManager.getInstance().getAminoAcids(Strand.POSITIVE, seq, defaultTable);
        assertEquals(3, acids.size());
        for (int i = 0; i < 3; i++) {
            assertEquals(aminoSeq[i], acids.get(i).getSymbol());
        }
    }

    @Test
    public void translateSequenceNegative() {
        // From URG$ "end" actually start of coding
        String seq = "ACGACGCCAT";   // "Extra" A at front (size == 10)
        char[] aminoSeq = {'S', 'A', 'M'};

        CodonTable defaultTable = CodonTableManager.getInstance().getDefaultCodonTable();
        List<CodonAA> acids = AminoAcidManager.getInstance().getAminoAcids(Strand.NEGATIVE, seq, defaultTable);
        assertEquals(3, acids.size());
        for (int i = 0; i < 3; i++) {
            assertEquals(aminoSeq[i], acids.get(i).getSymbol());
        }
    }

    @Test
    public void computeSequenceNegative() {
        // From URG$ "end" actually start of coding
        String seq = "ACGACGCCAT";   // "Extra" A at front (size == 10)
        char[] aminoSeq = {'S', 'A', 'M'};

        CodonTable defaultTable = CodonTableManager.getInstance().getDefaultCodonTable();
        AminoAcidSequence aaSequence = AminoAcidManager.getInstance().getAminoAcidSequence(Strand.NEGATIVE, 0, seq, defaultTable);
        assertEquals(3, aaSequence.getSequence().size());
        assertEquals(1, aaSequence.getStart());
        for (int i = 0; i < 3; i++) {
            assertEquals(aminoSeq[i], aaSequence.getSequence().get(i).getSymbol());
        }
    }

    //getAminoAcidSequence
    //GGCAGAACCAGCCGACGAGTCAGGCGCCGCATGGTCCCCTT
    @Test
    public void computeSequenceNegative2() {
        // From URG$ "end" actually start of coding
        String seq = "GGCAGAACCAGCCGACGAGTCAGGCGCCGCATGGTCCCCTT";
        byte[] aminoSeq1 = "LVLRRTLRRMTGK".getBytes();

        CodonTable defaultTable = CodonTableManager.getInstance().getDefaultCodonTable();
        AminoAcidSequence aaSequence = AminoAcidManager.getInstance().getAminoAcidSequence(Strand.NEGATIVE, 0, seq, defaultTable);
        assertEquals(13, aaSequence.getSequence().size());
        assertEquals(2, aaSequence.getStart());
        for (int i = 0; i < aminoSeq1.length; i++) {
            assertEquals(aminoSeq1[i], aaSequence.getSequence().get(i).getSymbol());
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
