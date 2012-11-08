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
        ResourceLocator locator = DBTable.createDBLocator(profile.getAbsolutePath());
        Connection conn = DBManager.getConnection(locator);
        Statement st = conn.createStatement();
        ResultSet rs = st.executeQuery("SELECT * FROM unigene");
        assertFalse(conn.isClosed());
        conn.close();
    }

}
