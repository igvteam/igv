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

        BasicFeature egfr = (BasicFeature) FeatureDB.getFeature(geneId);
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
        int exonNum = 2;

        BasicFeature lancl = (BasicFeature) FeatureDB.getFeature(geneId);
        Exon testExon = lancl.getExons().get(exonNum);
        Exon prevExon = lancl.getExons().get(exonNum-1);
        Exon nextExon = lancl.getExons().get(exonNum+1);
        AminoAcidSequence seq = testExon.getAminoAcidSequence(genome, prevExon, nextExon);
        assertEquals('I', seq.getSequence().get(1).getSymbol());

        AminoAcidManager.getInstance().setCodonTable(AminoAcidManager.DEFAULT_CODON_TABLE_PATH, 2);

        seq = testExon.getAminoAcidSequence(genome, prevExon, nextExon);
        assertEquals('H', seq.getSequence().get(1).getSymbol());

        AminoAcidSequence seq2 = testExon.getAminoAcidSequence(genome, prevExon, nextExon);

        //Shouldn't refetch
        assertEquals(seq, seq2);


    }

}
