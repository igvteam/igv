/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
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

    File f = new File(TestUtils.LARGE_DATA_DIR + "/snp130.bedz.sai");

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