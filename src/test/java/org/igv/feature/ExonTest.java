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

    /**
     * A leading codon split by the exon edge, with no neighboring exon to complete it, is skipped.  The rest of
     * the exon is translated in frame, the same as an exon that starts just past the split codon.
     */
    @Test
    public void testSplitCodonNoPrevExonPos() throws Exception {
        Exon cdsExon = ((BasicFeature) genome.getFeatureDB().getFeature("EGFR")).getExons().get(1);
        String chr = cdsExon.getChr();
        int start = cdsExon.getStart();
        int end = start + 30;

        Exon exon = new Exon(chr, start, end, Strand.POSITIVE);
        exon.setReadingFrame(2);    // 2 bases of the first codon are in the (missing) previous exon, 1 in this one
        Exon expected = new Exon(chr, start + 1, end, Strand.POSITIVE);
        expected.setReadingFrame(0);

        assertSameAminoAcids(expected.getAminoAcidSequence(genome, null, null), exon.getAminoAcidSequence(genome, null, null));
    }

    @Test
    public void testSplitCodonNoNextExonNeg() throws Exception {
        Exon cdsExon = ((BasicFeature) genome.getFeatureDB().getFeature("KRAS")).getExons().get(1);
        String chr = cdsExon.getChr();
        int start = cdsExon.getStart();
        int end = start + 30;

        Exon exon = new Exon(chr, start, end, Strand.NEGATIVE);
        exon.setReadingFrame(1);    // 1 base of the first codon is in the (missing) next exon, 2 in this one
        Exon expected = new Exon(chr, start, end - 2, Strand.NEGATIVE);
        expected.setReadingFrame(0);

        assertSameAminoAcids(expected.getAminoAcidSequence(genome, null, null), exon.getAminoAcidSequence(genome, null, null));
    }

    /**
     * A non-coding (UTR) neighbor cannot complete a split codon, it is treated the same as no neighbor.
     */
    @Test
    public void testSplitCodonNonCodingPrevExonPos() throws Exception {
        Exon cdsExon = ((BasicFeature) genome.getFeatureDB().getFeature("EGFR")).getExons().get(1);
        String chr = cdsExon.getChr();
        int start = cdsExon.getStart();
        int end = start + 30;

        Exon utr = new Exon(chr, start - 100, start - 50, Strand.POSITIVE);
        utr.setNonCoding(true);
        Exon exon = new Exon(chr, start, end, Strand.POSITIVE);
        exon.setReadingFrame(2);
        Exon expected = new Exon(chr, start + 1, end, Strand.POSITIVE);
        expected.setReadingFrame(0);

        assertSameAminoAcids(expected.getAminoAcidSequence(genome, null, null), exon.getAminoAcidSequence(genome, utr, null));
    }

    @Test
    public void testSplitCodonNonCodingNextExonNeg() throws Exception {
        Exon cdsExon = ((BasicFeature) genome.getFeatureDB().getFeature("KRAS")).getExons().get(1);
        String chr = cdsExon.getChr();
        int start = cdsExon.getStart();
        int end = start + 30;

        Exon utr = new Exon(chr, end + 50, end + 100, Strand.NEGATIVE);
        utr.setNonCoding(true);
        Exon exon = new Exon(chr, start, end, Strand.NEGATIVE);
        exon.setReadingFrame(1);
        Exon expected = new Exon(chr, start, end - 2, Strand.NEGATIVE);
        expected.setReadingFrame(0);

        assertSameAminoAcids(expected.getAminoAcidSequence(genome, null, null), exon.getAminoAcidSequence(genome, null, utr));
    }

    private void assertSameAminoAcids(AminoAcidSequence expected, AminoAcidSequence actual) {
        assertEquals(expected.getStart(), actual.getStart());
        assertEquals(expected.getSequenceString(), actual.getSequenceString());
    }

}
