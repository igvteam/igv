package org.igv.feature;

import org.igv.AbstractHeadlessTest;
import org.igv.feature.aa.AminoAcidSequence;
import org.igv.feature.aa.CodonTable;
import org.igv.feature.aa.CodonTableManager;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * User: jacob
 * Date: 2012/05/17
 */
public class ExonTest extends AbstractHeadlessTest {

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
     * @param geneId
     * @param exonNum
     * @param positive
     * @throws Exception
     */
    public void tstGetAminoAcidNumber(String geneId, int exonNum, boolean positive) throws Exception {

        int[] genomicOffsets = new int[]{0, 1, 2, 3, 4, 5};
        int[] expAANumbers = new int[]{1, 1, 1, 2, 2, 2};

        BasicFeature egfr = (BasicFeature) genome.getFeatureDB().getFeature(geneId);
        if (exonNum < 0) {
            exonNum = egfr.getExonCount() + exonNum;
        }
        Exon exon = egfr.getExons().get(exonNum);
        int ind = 0;
        int mult = positive ? 1 : -1;
        int start = positive ? exon.getCdStart() : exon.getCdEnd() - 1;
        for (Integer offset : genomicOffsets) {
            int genomicPosition = start + mult * offset;
            int AANumber = exon.getAminoAcidNumber(genomicPosition);
            assertEquals(expAANumbers[ind], AANumber);
            ind++;

        }

    }

    @Test
    public void testChangeCodonTable() throws Exception {

        String geneId = "LANCL2";
        int exonNum = 1;

        BasicFeature lancl = (BasicFeature) genome.getFeatureDB().getFeature(geneId);
        Exon testExon = lancl.getExons().get(exonNum);
        Exon prevExon = lancl.getExons().get(exonNum-1);
        Exon nextExon = lancl.getExons().get(exonNum+1);
        AminoAcidSequence seq = testExon.getAminoAcidSequence(genome, prevExon, nextExon);
        assertEquals('I', seq.getSequence().get(5).getSymbol());

        CodonTable codonTable= CodonTableManager.getInstance().getCodonTableByID(2);
        CodonTableManager.getInstance().setCurrentCodonTable(codonTable);
        seq = testExon.getAminoAcidSequence(genome, prevExon, nextExon);;
        assertEquals('M', seq.getSequence().get(5).getSymbol());

        AminoAcidSequence seq2 = testExon.getAminoAcidSequence(genome, prevExon, nextExon);

        //Shouldn't refetch
        assertEquals(seq, seq2);


    }

}
