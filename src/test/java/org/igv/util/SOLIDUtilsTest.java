package org.igv.util;

import static org.junit.Assert.assertEquals;
import org.junit.Test;
import org.igv.util.SOLIDUtils;

/**
 * User: jrobinso
 * Date: Feb 21, 2010
 */
public class SOLIDUtilsTest {
    @Test
    public void testConvertToColorSpace() {
        // Example from Solid data mgt app node


        byte[] sequence = "GTGCACCGTGCACG".getBytes();
        byte[] expCS = new byte[]{'G', 1, 1, 3, 1, 1, 0, 3, 1, 1, 3, 1, 1, 3};

        byte[] cs = SOLIDUtils.convertToColorSpace(sequence);
        for (int i = 0; i < sequence.length; i++) {
            assertEquals(String.valueOf(i), expCS[i], cs[i]);
        }
        // Add your code here
    }
}
