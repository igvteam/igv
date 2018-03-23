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

package org.broad.igv.maf;

import org.broad.igv.util.TestUtils;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.util.List;

import static junit.framework.Assert.assertEquals;

/**
 * @author jrobinso
 *         Date: 2/18/13
 *         Time: 11:48 AM
 */
public class MAFParserTest {
    @Test
    public void testParse() throws Exception {

        String mafFile = TestUtils.DATA_DIR + "maf/ucscSample.maf";

        MAFIndex.blockSize = 1;

        MAFParser parser = new MAFParser(mafFile);

        List<MultipleAlignmentBlock> alignments = parser.loadAlignments("chr1", 0, 1000000);

        assertEquals(3, alignments.size());

        // Last alignment has 2 gaps of 3 bases in human relative ot guinea pig
        // s hg18.chr1                       43219 50 + 247249719 G---ATG---TCATAATAAATGGTGCATATCCAGAGTGCAAGATGATTCAGTCTCA
        // s cavPor3.scaffold_32          16022892 55 +  25641242 GGACATGGAACAATAATAACTGATGCAC-TGTAGAGCACAATATGATATAGTTTCT

        MultipleAlignmentBlock ma = alignments.get(2);
        assertEquals(43219, ma.getStart());
        assertEquals(43219 + 50, ma.getEnd());
        List<MultipleAlignmentBlock.Gap> gaps = ma.getGaps();
        assertEquals(2, gaps.size());

        MultipleAlignmentBlock.Gap firstGap = gaps.get(0);
        assertEquals(43219 + 1, firstGap.position, 1.0e-6);
        assertEquals(1, firstGap.startIdx);
        assertEquals(3, firstGap.size);

        MultipleAlignmentBlock.Gap secondGap = gaps.get(1);
        assertEquals(43219 + 4, secondGap.position, 1.0e-6);
        assertEquals(7, secondGap.startIdx);
        assertEquals(3, secondGap.size);

        for (int i = ma.getStart(); i < ma.getEnd(); i++) {
            System.out.println(i + "  " + ma.getGapAdjustedIndex(i));
        }

        assertEquals("hg18 Multiz", parser.getTrackName());

        String[] expectedSpeciesOrder = new String[]{"hg18", "tarSyr1", "gorGor1", "panTro2", "ponAbe2", "bosTau4",
                "equCab2", "ornAna1", "cavPor3", "canFam2", "rheMac2"};
        Object [] species = parser.getSpecies().toArray();
        Assert.assertArrayEquals(expectedSpeciesOrder, species);

        String indexFile = mafFile + ".index";
        (new File(indexFile)).delete();

    }
}
