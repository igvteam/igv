package org.broad.igv.feature.tribble;

import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.Feature;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * @author jrobinso
 *         Date: 1/23/14
 *         Time: 4:24 PM
 */
public class EncodePeakCodecTest {

    @Test
    public void testDecode() throws Exception {

        String line = "chr22\t16847690\t16857344\t.\t354\t.\t5.141275\t14.1\t-1\t100";
        EncodePeakCodec codec = new EncodePeakCodec();

        AsciiFeatureCodec codec2 = new EncodePeakCodec();

        assertTrue(codec.getClass() == codec2.getClass());

        Feature feature =  codec.decode(line);
        Feature feature2 = codec2.decode(line);

        assertEquals("chr22", feature.getChr());


    }


}
