package org.igv.feature.tribble;

import org.igv.AbstractHeadlessTest;
import org.igv.feature.IGVFeature;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * User: jacob
 * Date: 2012-Dec-11
 */

public class IGVBEDCodecTest extends AbstractHeadlessTest {

    @Test
    public void testTabDelimited() throws Exception {
        String line = "chr1\t1051161\t1051177\tstrong GGGCGGGTGGGGCGGG\t0.81";
        IGVBEDCodec codec = new IGVBEDCodec();
        IGVFeature feature =  codec.decode(line);
        assertEquals("strong GGGCGGGTGGGGCGGG", feature.getName());
    }

    @Test
    public void testMultiTabDelimited() throws Exception {
        String line = "chr1\t1051161\t\t\t1051177\tstrong GGGCGGGTGGGGCGGG\t0.81";
        IGVBEDCodec codec = new IGVBEDCodec();
        IGVFeature feature =  codec.decode(line);
        assertEquals("strong GGGCGGGTGGGGCGGG", feature.getName());
    }

    @Test
    public void testSpaceDelimited() throws Exception {
        String line = "chr1 1051161 1051177 strong_GGGCGGGTGGGGCGGG 0.81";
        IGVBEDCodec codec = new IGVBEDCodec();
        IGVFeature feature =  codec.decode(line);
        assertEquals("strong_GGGCGGGTGGGGCGGG", feature.getName());
    }

}
