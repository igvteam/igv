package org.broad.igv.feature.exome;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.ui.panel.ReferenceFrame;
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
public class ExomeUtilsTest extends AbstractHeadlessTest {

    static ExomeReferenceFrame frame;

    @BeforeClass
    public static void setup() throws IOException {
        Map<String, List<Feature>> allFeatures = loadTestFeatures();

        frame = new ExomeReferenceFrame(new ReferenceFrame("test"), allFeatures);
        frame.jumpTo(new Locus("chr6", 1, 100));
    }



    @Test
    public void genomeToExomePosition() {
        final int genomeStart = 10831323;
        final int exomeStart = 2685;

        // Genome positino at exact start of blocks
        int calcExomePosition = frame.genomeToExomePosition(genomeStart);
        assertEquals(exomeStart, calcExomePosition);

        // Genome position in interior of blocks
        int genomePosition = genomeStart + 100;
        int expectedExomePosition = exomeStart + 100;
        calcExomePosition = frame.genomeToExomePosition(genomePosition);
        assertEquals(expectedExomePosition, calcExomePosition);

        // Between 2 blocks -- position exome at end of first blocks
        //Block 2 [2208833, 2208917, 223, 84]
        //Block 3 [2214930, 2215032, 307, 102]
        genomePosition = (2208917 + 2214930) / 2;
        expectedExomePosition = 223 + 84;
        calcExomePosition = frame.genomeToExomePosition(genomePosition);
        assertEquals(expectedExomePosition, calcExomePosition);

        // Before first blocks
        genomePosition = 100;
        expectedExomePosition = 0;
        calcExomePosition = frame.genomeToExomePosition(genomePosition);
        assertEquals(expectedExomePosition, calcExomePosition);

        // After last blocks
        //Block 100 [111693771, 111696954, 29865, 3183]
        genomePosition = 111696954 + 100;
        expectedExomePosition = 29865 + 3183;
        calcExomePosition = frame.genomeToExomePosition(genomePosition);
        assertEquals(expectedExomePosition, calcExomePosition);
    }


    @Test
    public void exomeToGenomePosition() {
        final int genomeStart = 10831323;
        final int exomeStart = 2685;

        int exomePosition = exomeStart + 100;
        int expectedGenomePosition = genomeStart + 100;

        int calcGenomePosition = frame.exomeToGenomePosition(exomePosition);
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
