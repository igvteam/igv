package org.broad.igv.feature.tribble;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @author jrobinso
 *         Date: 4/8/13
 *         Time: 2:51 PM
 */
public class MUTCodecTest extends AbstractHeadlessTest {

    String path = TestUtils.DATA_DIR + "maf/TCGA_GBM_Level3_Somatic_Mutations_08.28.2008.maf";



    @Test
    public void testIsMutationAnnotationFile() throws Exception {
        assertTrue(MUTCodec.isMutationAnnotationFile(path));
    }

    @Test
    public void testGetSamples() throws Exception {

        int expectedSampleCount = 146;
        MUTCodec codec = new MUTCodec(path, genome);
        String [] samples = codec.getSamples();
        assertEquals(expectedSampleCount, samples.length);
        assertTrue(samples[0].startsWith("TCGA"));

    }


}
