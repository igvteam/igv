package org.broad.igv.sam;

import org.junit.Test;
import static org.junit.Assert.*;

public class ByteSubarrayTest {

    @Test
    public void testReverseComplement()  {

        byte [] forwardSeq = "AAAAATCGAAAA".getBytes();
        byte [] expected = "CGAT".getBytes();

        ByteSubarray subarray = new ByteSubarray(forwardSeq, 4, 4, (byte) 0);
        ByteSubarray reverse = AlignmentUtils.reverseComplementCopy(subarray);
        for(int i=0; i<4; i++) {
            byte b = reverse.getByte(i);
            assertEquals(expected[i], b);
        }
    }
}