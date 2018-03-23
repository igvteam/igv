/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
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