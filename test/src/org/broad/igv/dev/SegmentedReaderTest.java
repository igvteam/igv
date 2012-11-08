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

package org.broad.igv.dev;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.data.seg.SegmentedAsciiDataSet;
import org.broad.igv.dev.db.DBManager;
import org.broad.igv.dev.db.DBTable;
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
        DBTable table = DBTable.build(locator, tableName);
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
        DBTable table = DBTable.build(locator, tableName);

        SegmentedReader reader = new SegmentedReader(locator);
        SegmentedAsciiDataSet ds = reader.loadFromDB(table);

    }

}
