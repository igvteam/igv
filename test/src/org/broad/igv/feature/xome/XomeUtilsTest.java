package org.broad.igv.feature.xome;

import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
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
public class XomeUtilsTest {

    @Test
    public void testCollapseTranscripts() throws Exception {

        Map<String, List<Feature>> allFeatures = loadTestFeatures();

        List<Feature> chr1Features = allFeatures.get("chr1");

        List<Block> blocks = XomeUtils.collapseTranscripts(chr1Features);

        int lastEnd = -1;
        int exomeStart = 0;
        int blockIdx = 0;
        for(Block b : blocks) {
            assertTrue(b.getGenomeEnd() > lastEnd);
            assertEquals(b.getExomeStart(), exomeStart);
            assertEquals(blockIdx, b.getIdx());

            lastEnd = b.getGenomeEnd();
            exomeStart += b.getLength();
            blockIdx++;
        }

    }


    Map<String, List<Feature>> loadTestFeatures() throws IOException {

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
