package org.broad.igv.feature.exome;

import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @author Jim Robinson
 * @date 5/24/12
 */
public class ExomeUtilsTest {

    static List<ExomeBlock> chr6Blocks;

    @BeforeClass
    public static void setup() throws IOException {
        Map<String, List<Feature>> allFeatures = loadTestFeatures();
        chr6Blocks = ExomeUtils.collapseTranscripts(allFeatures.get("chr6"));

    }

    @Test
    public void testCollapseTranscripts() throws Exception {
        int lastGenomeEnd = 0;
        int lastExomeEnd = 0;
        int exomeStart = 0;
        int blockIdx = 0;
        for (ExomeBlock b : chr6Blocks) {
            assertTrue(b.getGenomeEnd() > lastGenomeEnd);
            assertEquals(exomeStart, b.getExomeStart());
            assertEquals(exomeStart, lastExomeEnd);
            assertEquals(blockIdx, b.getIdx());

            lastGenomeEnd = b.getGenomeEnd();
            lastExomeEnd = b.getExomeEnd();
            exomeStart += b.getLength();
            blockIdx++;
        }

    }

    //uc010joq.1	chr6	+	10831323	10837542	10832832	10837061	5	10831323,10832788,10833179,10834125,10836858,	10831460,10832852,10833256,10834227,10837542,	Q5T4I6	uc010joq.1
//    @Test
//    public void testGetByGenomePosition() {
//
//        final int genomePosition = 10831324;
//        ExomeBlock b = ExomeUtils.getBlockAtGenomePosition(chr6Blocks, genomePosition);
//        assertTrue(genomePosition >= b.getGenomeStart() && genomePosition < b.getGenomeEnd());
//    }
//
//    @Test
//    public void testGetByExomePosition() {
//
//        final int exomePosition = 3000;
//        ExomeBlock b = ExomeUtils.getBlockAtExomePosition(chr6Blocks, exomePosition);
//        assertTrue(exomePosition >= b.getExomeStart() && exomePosition < b.getExomeEnd());
//    }

    @Test
    public void genomeToExomePosition() {
        final int genomeStart = 10831323;
        final int exomeStart = 2685;

        // Genome positino at exact start of block
        int calcExomePosition = ExomeUtils.genomeToExomePosition(chr6Blocks, genomeStart);
        assertEquals(exomeStart, calcExomePosition);

        // Genome position in interior of block
        int genomePosition = genomeStart + 100;
        int expectedExomePosition = exomeStart + 100;
        calcExomePosition = ExomeUtils.genomeToExomePosition(chr6Blocks, genomePosition);
        assertEquals(expectedExomePosition, calcExomePosition);

        // Between 2 blocks -- position exome at end of first block
        //Block 2 [2208833, 2208917, 223, 84]
        //Block 3 [2214930, 2215032, 307, 102]
        genomePosition = (2208917 + 2214930) / 2;
        expectedExomePosition = 223 + 84;
        calcExomePosition = ExomeUtils.genomeToExomePosition(chr6Blocks, genomePosition);
        assertEquals(expectedExomePosition, calcExomePosition);

        // Before first block
        genomePosition = 100;
        expectedExomePosition = 0;
        calcExomePosition = ExomeUtils.genomeToExomePosition(chr6Blocks, genomePosition);
        assertEquals(expectedExomePosition, calcExomePosition);

        // After last block
        //Block 100 [111693771, 111696954, 29865, 3183]
        genomePosition = 111696954 + 100;
        expectedExomePosition = 29865 + 3183;
        calcExomePosition = ExomeUtils.genomeToExomePosition(chr6Blocks, genomePosition);
        assertEquals(expectedExomePosition, calcExomePosition);
    }

    @Test
    public void testGetBlockByGenomePosition() {
        int genomePosition = (2208917 + 2214930) / 2;
        int idx = ExomeUtils.getIndexForGenomePosition(chr6Blocks, genomePosition);
        assertEquals(2, idx);

        ExomeBlock b = chr6Blocks.get(idx);
        assertTrue(genomePosition >= b.getGenomeStart());

        ExomeBlock nextBlock = chr6Blocks.get(idx + 1);
        assertTrue(genomePosition < nextBlock.getGenomeStart());

    }

    @Test
    public void exomeToGenomePosition() {
        final int genomeStart = 10831323;
        final int exomeStart = 2685;

        int exomePosition = exomeStart + 100;
        int expectedGenomePosition = genomeStart + 100;

        int calcGenomePosition = ExomeUtils.exomeToGenomePosition(chr6Blocks, exomePosition);
        assertEquals(expectedGenomePosition, calcGenomePosition);
    }


    static Map<String, List<Feature>> loadTestFeatures() throws IOException {

        Map<String, List<Feature>> allFeatures = new HashMap<String, List<Feature>>();

        String file = TestUtils.DATA_DIR + "gene/UCSCgenes_sample.gene";
        FeatureCodec codec = CodecFactory.getCodec(file, null);
        AbstractFeatureReader<Feature> bfs = AbstractFeatureReader.getFeatureReader(file, codec, false);
        Iterable<Feature> iter = bfs.iterator();
        for (Feature f : iter) {
            List<Feature> flist = allFeatures.get(f.getChr());
            if (flist == null) {
                flist = new ArrayList<Feature>(5000);
                allFeatures.put(f.getChr(), flist);
            }
            flist.add(f);
        }

        return allFeatures;

    }
}
