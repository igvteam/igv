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

package org.broad.igv.dev.db;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.After;
import org.junit.Test;

import java.io.File;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.Statement;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertFalse;


/**
 * User: jacob
 * Date: 2012/05/29
 */
public class DBManagerTest extends AbstractHeadlessTest {

    @After
    public void tearDown() throws Exception {
        super.tearDown();
        DBManager.shutdown();
    }

    static ResultSet getAllFromSQLTable(String path, String table) throws Exception {
        File dataDir = new File(TestUtils.DATA_DIR);
        String url = DBManager.createConnectionURL("sqlite", dataDir.getAbsolutePath(), path, null);
        ResourceLocator locator = new ResourceLocator(url);
        Connection conn = DBManager.getConnection(locator);

        String query = "SELECT * FROM " + table;
        Statement st = conn.createStatement();
        return st.executeQuery(query);
    }

    @Test
    public void testConnectSQLite() throws Exception {
        String path = "sql/unigene.db";
        ResultSet rs = getAllFromSQLTable(path, "unigene");

        assertEquals(12, rs.getMetaData().getColumnCount());
        rs.close();
    }

    @Test
    public void testGetStoredConnection() throws Exception {
        String subPath = "sql/unigene_profile.dbxml";
        File profile = new File(TestUtils.DATA_DIR, subPath);
        ResourceLocator locator = DBProfile.parseProfile(profile.getAbsolutePath()).getDBLocator();
        Connection conn = DBManager.getConnection(locator);
        Statement st = conn.createStatement();
        ResultSet rs = st.executeQuery("SELECT * FROM unigene");
        assertFalse(conn.isClosed());
        conn.close();
    }

}
