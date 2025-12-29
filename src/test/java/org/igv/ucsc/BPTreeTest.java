package org.igv.ucsc;

import org.igv.util.TestUtils;
import org.junit.Test;

import java.io.IOException;
import java.nio.ByteOrder;

import static org.junit.Assert.*;

public class BPTreeTest {

    @Test
    public void testSearch() throws IOException {
        String testFile = TestUtils.DATA_DIR + "twobit/GCA_004363605.1.2bit.bpt";

        BPTree tree =  BPTree.loadBPTree(testFile, 0);

        assertNotNull(tree);
        assertEquals(256, tree.blockSize);
        assertEquals(15, tree.keySize);
        assertEquals(8, tree.valSize);

        long result = tree.searchForOffset("RJWJ011179649.1");
        assertEquals(748759988, result);
        assertNotNull(result);
    }

//    @Test
//    public void testSearch2() throws IOException {
//        String testFile = TestUtils.DATA_DIR + "twobit/GCA_004363605.1.2bit";
//
//        SeekableStream is = IGVSeekableStreamFactory.getInstance().getStreamFor(testFile);
//
//        TwoBitIndex tree = new TwoBitIndex(is, ByteOrder.LITTLE_ENDIAN, );
//
//        assertNotNull(tree);
//        assertEquals(256, tree.blockSize);
//        assertEquals(15, tree.keySize);
//        assertEquals(8, tree.valSize);
//
//        long[] result = tree.search("RJWJ011179649.1");
//        assertEquals(748759988, result[0]);
//        assertNotNull(result);
//    }


}