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
 * Class for loading ALL of the data from a single table
 * in a single database.
 *
 * @author Jacob Silterra
 * @date 2012/05/30
 */
public abstract class WholeTableDBReader<T> {

    private static Logger log = Logger.getLogger(WholeTableDBReader.class);

    protected ResourceLocator locator;
    protected String tableName;
    protected String queryString = "SELECT * FROM ";

    public WholeTableDBReader(ResourceLocator locator, String tableName) {
        this.locator = locator;
        this.tableName = tableName;
        queryString += tableName;
    }

    public T load() {

        Connection conn = null;
        Statement st = null;
        ResultSet rs = null;
        T obj = null;

        try {
            conn = DBManager.getConnection(locator);
            st = conn.createStatement();
            rs = st.executeQuery(queryString);

            obj = processResultSet(rs);

        } catch (SQLException e) {
            log.error("Database error", e);
            throw new RuntimeException("Database error", e);
        } finally {
            DBManager.closeResources(rs, st, conn);
        }

        return obj;
    }


    protected abstract T processResultSet(ResultSet rs) throws SQLException;

}
