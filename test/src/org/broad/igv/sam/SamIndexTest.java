/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.sam;

import org.broad.igv.sam.reader.FeatureIndex;
import org.broad.igv.util.TestUtils;
import org.junit.AfterClass;
import org.junit.BeforeClass;
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
    @Test
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