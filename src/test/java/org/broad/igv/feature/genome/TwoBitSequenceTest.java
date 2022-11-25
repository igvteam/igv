package org.broad.igv.feature.genome;

import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.IOException;

import static org.junit.Assert.*;

public class TwoBitSequenceTest {



    @Test
    public void getSequence() throws IOException {
        String testFile = TestUtils.DATA_DIR + "twobit/foo.2bit";
        TwoBitSequence seq = new TwoBitSequence(testFile);
        assertNotNull(seq);

    }
}