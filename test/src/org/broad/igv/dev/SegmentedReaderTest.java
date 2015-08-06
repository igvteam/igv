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

package org.broad.igv.dev;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.data.seg.SegmentedAsciiDataSet;
import org.broad.igv.dev.db.DBManager;
import org.broad.igv.dev.db.DBProfile;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.File;

import static junit.framework.Assert.*;

/**
 * Test of our multi-purpose reader class
 * User: jacob
 * Date: 2012-Oct-23
 */
public class SegmentedReaderTest extends AbstractHeadlessTest {


    private SegmentedAsciiDataSet loadSegFile() throws Exception {
        String filePath = TestUtils.DATA_DIR + "seg/canFam2_hg18.seg";
        ResourceLocator locator = new ResourceLocator(filePath);
        SegmentedReader reader = new SegmentedReader(locator);
        SegmentedAsciiDataSet ds = reader.loadFromFile();

        assertEquals(18, ds.getChromosomes().size());

        return ds;
    }

    @Test
    public void testLoadSegFile() throws Exception {
        SegmentedAsciiDataSet ds = loadSegFile();
        assertNotNull(ds);
        assertFalse(ds.getChromosomes().isEmpty());
    }

    @Test
    public void testLoadSegDB() throws Exception {
        SegmentedAsciiDataSet ds = loadSegDB();
        assertNotNull(ds);
        assertFalse(ds.getChromosomes().isEmpty());
    }

    private SegmentedAsciiDataSet loadSegDB() throws Exception {
        String host = (new File(TestUtils.DATA_DIR)).getAbsolutePath();
        String path = "seg/canFam2_hg18.db";
        String url = DBManager.createConnectionURL("sqlite", host, path, null);
        ResourceLocator locator = new ResourceLocator(url);
        String tableName = "canFam2";

        SegmentedReader reader = new SegmentedReader(locator);
        DBProfile.DBTable table = DBProfile.DBTable.build(locator, tableName);
        SegmentedAsciiDataSet ds = reader.loadFromDB(table);

        return ds;
    }

    @Test
    public void compareFileDB() throws Exception {
        SegmentedAsciiDataSet fileDS = loadSegFile();
        SegmentedAsciiDataSet dbDS = loadSegDB();

        assertEquals(fileDS.getChromosomes(), dbDS.getChromosomes());
        assertEquals(fileDS.getSampleNames(), dbDS.getSampleNames());
    }

    //Scratch work for testing loading from local mysql db
    public void testLoadMySQL() throws Exception {
        String url = DBManager.createConnectionURL("mysql", "calcium", "igv_nobel_dev", null);
        ResourceLocator locator = new ResourceLocator(url);
        locator.setUsername("igv_nobel_dev");
        locator.setPassword("nottherealpassword");
        String tableName = "CNV";
        DBProfile.DBTable table = DBProfile.DBTable.build(locator, tableName);

        SegmentedReader reader = new SegmentedReader(locator);
        SegmentedAsciiDataSet ds = reader.loadFromDB(table);

    }

}
