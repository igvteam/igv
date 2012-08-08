package org.broad.igv.hic.matrix;

import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * @author jrobinso
 *         Date: 8/8/12
 *         Time: 8:59 AM
 */
public class DiskResidentBlockMatrixTest {

    static DiskResidentBlockMatrix blockMatrix;

    @BeforeClass
    public static void setUp() throws Exception {

        blockMatrix = new DiskResidentBlockMatrix("/Users/jrobinso/projects/hic/block_chr14.bin");

    }

    @Test
    public void testInit() throws Exception {

        assertEquals("14", blockMatrix.getChr1());

        assertEquals(8, blockMatrix.getRemSize());

    }

    @Test
    public void testLoadLastBlock() throws Exception {

        int rowBlockIdx = 10;
        int rowColIdx = 7;

        float [][] data = blockMatrix.loadBlockData(rowBlockIdx, rowColIdx);

        assertEquals(8, data.length);
        assertEquals(10, data[0].length);


    }
}
