package org.broad.igv.feature.genome;


import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by jrobinso on 6/13/17.
 */
public class TwoBitSequenceTest extends AbstractHeadlessTest {

    String path = "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit" ; //TestUtils.DATA_DIR + "fasta/eboVir3.2bit";

    @Test
    public void testInit() throws Exception {

        TwoBitSequence seq = new TwoBitSequence(path);

        assertTrue(seq  != null);
    }


//    @Test
//    public void getSequence() throws Exception {
//
//    }
//
//    @Test
//    public void getBase() throws Exception {
//
//    }
//
//    @Test
//    public void getChromosomeNames() throws Exception {
//
//    }
//
//    @Test
//    public void getChromosomeLength() throws Exception {
//
//    }

}
