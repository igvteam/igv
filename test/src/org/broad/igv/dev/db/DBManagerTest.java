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
import org.broad.igv.util.TestUtils;
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


    @Test
    public void testConnectSQLite() throws Exception{
        String path = "sql/unigene.db";
        File dataDir = new File(TestUtils.DATA_DIR);
        String url = DBManager.createConnectionURL(dataDir.getAbsolutePath(), path, null, "sqlite");
        Connection conn = DBManager.getConnection(url);

        String query = "SELECT * FROM unigene";
        Statement st = conn.createStatement();
        ResultSet rs = st.executeQuery(query);

        int count = 0;
        assertEquals(12, rs.getMetaData().getColumnCount());

    }

    @Test
    public void testGetStoredConnection() throws Exception{
        String subPath = "sql/unigene_profile.xml";
        File profile = new File(TestUtils.DATA_DIR, subPath);
        Connection conn = DBManager.getStoredConnection(profile.getAbsolutePath());
        assertFalse(conn.isClosed());
    }

}
