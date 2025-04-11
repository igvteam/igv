package org.broad.igv.ucsc.bb;

import org.broad.igv.ucsc.BPTree;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.IOException;
import java.util.Map;

import static org.junit.Assert.*;

public class ChromTreeTest {

    @Test
    public void testBPTree() throws IOException {

        String bbFile = TestUtils.DATA_DIR + "bb/cytoBandMapped.bb";

        // BBFile file = new BBFile(bbFile, null);
        ChromTree chromTree = new ChromTree(bbFile, 757);
        assertNotNull(chromTree);

        int chrID = chromTree.getIdForName("chr10");
        String chrName = chromTree.getNameForId(chrID);
        assertEquals("chr10", chrName);
    }

    /**
     * Test querying a bb file with a very large chrom tree.
     */
    @Test
    public void testLargeTree() throws IOException {

        String bbFile = "https://hgdownload.soe.ucsc.edu/hubs/GCA/004/027/955/GCA_004027955.1/GCA_004027955.1.chromAlias.bb";

        // BBFile file = new BBFile(bbFile, null);
        ChromTree chromTree = new ChromTree(bbFile, 738);
        assertNotNull(chromTree);
        assertEquals(7677333, chromTree.getItemCount());
        int chrID = chromTree.getIdForName("PVKH010020857.1");
        String chrName = chromTree.getNameForId(chrID);
        assertEquals("PVKH010020857.1", chrName);
    }

    @Test
    public void testChromSizeEstimate() throws IOException {
        String bbFile = "https://hgdownload.soe.ucsc.edu/gbdb/hs1/ncbiRefSeq/ncbiRefSeqCurated.bb";
        ChromTree chromTree = new ChromTree(bbFile, 1752);
        assertNotNull(chromTree);
        long estSize = chromTree.estimateGenomeSize();
        assertEquals(3117275501l, estSize);
    }

    /**
     * Test a BB file with a very large chrom tree (> 7 million contigs
     */
    @Test
    public void testChromSizeEstimate2() throws IOException {

        String bbFile = "https://hgdownload.soe.ucsc.edu/hubs/GCA/004/027/955/GCA_004027955.1/GCA_004027955.1.chromAlias.bb";
        ChromTree chromTree = new ChromTree(bbFile, 738);
        System.out.println(chromTree.getItemCount());
        assertNotNull(chromTree);
        long estSize = chromTree.estimateGenomeSize();
        assertTrue(estSize > 0);

    }
}