/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.igv.sam;

import org.igv.sam.reader.FeatureIndex;
import org.igv.util.TestUtils;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import java.io.File;

import static junit.framework.Assert.assertTrue;

/**
 * @author jrobinso
 */
public class SamIndexTest {

    File f = new File(TestUtils.LARGE_DATA_DIR + "snp130.bedz.sai");

    public SamIndexTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    /**
     * Test of getIndexedChromosomes method, of class SamIndex.
     */
    @Test @Ignore("Requires largedata bundle")
    public void testGetIndexedChromosomes() {
        FeatureIndex idx = new FeatureIndex(f);
        assertTrue(idx.getIndexedChromosomes().size() > 0);
    }

    /**
     * Test of add method, of class SamIndex.
     */
    @Test
    public void testAdd() {
    }

    /**
     * Test of getTileDef method, of class SamIndex.
     */
    @Test
    public void testGetTileDef() {
    }

    /**
     * Test of store method, of class SamIndex.
     */
    @Test
    public void testStore() throws Exception {
    }

    /**
     * Test of getTileWidth method, of class SamIndex.
     */
    @Test
    public void testGetTileWidth() {
    }

}