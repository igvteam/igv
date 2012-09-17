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

        List<AminoAcid> acids = AminoAcidManager.getInstance().getAminoAcids(seq, Strand.POSITIVE);
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

        List<AminoAcid> acids = AminoAcidManager.getInstance().getAminoAcids(seq, Strand.NEGATIVE);
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
            assertTrue(AminoAcidManager.getInstance().setCodonTableById(id));
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
            boolean loaded = aam.setCodonTableById(codonTableIds[ii]);

            assertTrue("Failed to load codon table with id " + codonTableIds[ii], loaded);

            AminoAcid actualAA = aam.getAminoAcid(testCodons[ii]);
            assertNotSame("Got null amino acid for " + testCodons[ii], AminoAcid.NULL_AMINO_ACID, actualAA);
            assertEquals(String.valueOf(expAAs[ii]), String.valueOf(actualAA.getSymbol()));

            //We want to only store one copy of each amino acid. Check that this is the case
            assertTrue(AminoAcidManager.getAminoAcidByName(actualAA.getShortName()) == actualAA);
        }
    }

//    public static void main(String[] args) throws Exception{
//        //Test loading asn.1 using bouncycastle
    //binary asn.1 data from ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.val
    //bouncycastle from http://www.bouncycastle.org/
//        String filePath = TestUtils.DATA_DIR + "/gc.val";
//        InputStream is = new FileInputStream(filePath);
//        ASN1InputStream ASNis = new ASN1InputStream(is);
//        ASN1Primitive obj = ASNis.readObject();
//        BERSet set = (BERSet) obj;
//        ASN1Encodable[] objs = set.toArray();
//        for(ASN1Encodable aobj: objs){
//            byte[] data = aobj.toASN1Primitive().getEncoded();
//            System.out.println(new String(data));
//        }
//    }
}
