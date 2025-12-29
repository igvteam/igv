package org.broad.igv.sam.smrt;

import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class SMRTKineticsDecoderTest {

    public SMRTKineticsDecoderTest() {
    }

    @Test
    public void testDecoder() throws Exception {
 
        short v0 = SMRTKineticsDecoder.lookupFrameCount((byte) 0);
        short v128 = SMRTKineticsDecoder.lookupFrameCount((byte) 128);
        short v255 = SMRTKineticsDecoder.lookupFrameCount((byte) 255);

        assertEquals(0, v0);
        assertEquals(192, v128);
        assertEquals(952, v255);
    }
}