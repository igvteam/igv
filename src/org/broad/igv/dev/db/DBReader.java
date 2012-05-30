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

import org.apache.log4j.Logger;
import org.broad.igv.util.ResourceLocator;

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

/**
 * @author Jim Robinson
 * @date 1/20/12
 */
public abstract class DBReader<T> {

    private static Logger log = Logger.getLogger(DBReader.class);


    public static void closeResources(ResultSet rs, Statement st, Connection conn) {
        if (rs != null) {
            try {
                rs.close();
            } catch (SQLException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }
        if (st != null) {
            try {
                st.close();
            } catch (SQLException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }
        if (conn != null) try {
            conn.close();
        } catch (SQLException e) {
            log.error("Error closing sql connection", e);
        }

    }

    public T load(ResourceLocator locator, String sql) {


        ResultSet rs = null;
        Statement st = null;
        Connection conn = null;
        T obj = null;
        try {
            conn = DBManager.getConnection(locator);
            st = conn.createStatement();
            rs = st.executeQuery(sql);

            obj = processResultSet(rs);

            System.out.println("Disconnected from database");
        } catch (SQLException e) {
            log.error("Database error", e);
            throw new RuntimeException("Database error", e);
        } finally {
            closeResources(rs, st, conn);
        }

        return obj;
    }

    protected abstract T processResultSet(ResultSet rs) throws SQLException;
}
