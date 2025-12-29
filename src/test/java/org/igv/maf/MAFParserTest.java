package org.igv.maf;

import org.igv.util.TestUtils;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.util.List;

import static junit.framework.Assert.assertEquals;

/**
 * @author jrobinso
 * Date: 2/18/13
 * Time: 11:48 AM
 */
public class MAFParserTest {
    @Test
    public void testParse() throws Exception {

        String mafFile = TestUtils.DATA_DIR + "maf_ucsc/ucscSample.maf";

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
        assertEquals("hg18 Multiz", parser.getTrackName());

        String[] expectedSpeciesOrder = new String[]{"hg18", "tarSyr1", "gorGor1", "panTro2", "ponAbe2", "bosTau4",
                "equCab2", "ornAna1", "cavPor3", "canFam2", "rheMac2"};
        Object[] species = parser.getSpecies().toArray();
        Assert.assertArrayEquals(expectedSpeciesOrder, species);

        String indexFile = mafFile + ".index";
        (new File(indexFile)).delete();

    }
}
