package org.broad.igv.maf;

import com.mongodb.util.MyAsserts;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.util.List;

import static com.mongodb.util.MyAsserts.assertTrue;
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

        MAFParser parser = new MAFParser(mafFile);

        List<MultipleAlignmentBlock> alignments = parser.parse();

        assertTrue(alignments.size() == 3);

        // Last alignment has 2 gaps of 3 bases in human relative ot guinea pig
        // s hg18.chr1                       43219 50 + 247249719 G---ATG---TCATAATAAATGGTGCATATCCAGAGTGCAAGATGATTCAGTCTCA
        // s cavPor3.scaffold_32          16022892 55 +  25641242 GGACATGGAACAATAATAACTGATGCAC-TGTAGAGCACAATATGATATAGTTTCT

        MultipleAlignmentBlock ma = alignments.get(2);
        assertEquals(43219, ma.getStart());
        assertEquals(43219 + 50, ma.getEnd());
        List<MultipleAlignmentBlock.Gap> gaps = ma.getGaps();
        assertEquals(2, gaps.size());

        MultipleAlignmentBlock.Gap firstGap = gaps.get(0);
        assertEquals(43219 + 0.5, firstGap.position, 1.0e-6);
        assertEquals(1, firstGap.startIdx);
        assertEquals(3, firstGap.size);

        MultipleAlignmentBlock.Gap secondGap = gaps.get(1);
        assertEquals(43219 + 3.5, secondGap.position, 1.0e-6);
        assertEquals(7, secondGap.startIdx);
        assertEquals(3, secondGap.size);

        for(int i=ma.getStart(); i<ma.getEnd(); i++) {
            System.out.println(i + "  " + ma.getGapAdjustedIndex(i));
        }



    }
}
