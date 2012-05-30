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
import org.broad.tribble.CloseableTribbleIterator;

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Iterator;

/**
 * Class for reading only portions of a table (queries) repeatedly.
 * The connection is NOT closed between queries.
 *
 * @author Jacob Silterra
 * @date 2012/05/30
 */
public abstract class DBReader<T> {

    private static Logger log = Logger.getLogger(DBReader.class);

    protected ResourceLocator locator;
    protected String table;
    protected String baseQueryString = "SELECT * FROM ";

    public DBReader(ResourceLocator locator, String table) {
        this.locator = locator;
        this.table = table;
        baseQueryString += table;
    }

    protected ResultSet loadResultSet(String queryString) {

        try {
            Connection conn = DBManager.getConnection(locator);
            Statement st = conn.createStatement();
            ResultSet rs = st.executeQuery(queryString);

            return rs;

        } catch (SQLException e) {
            log.error("Database error", e);
            throw new RuntimeException("Database error", e);
        }
    }

    protected CloseableTribbleIterator loadIterator(String queryString) {
        return new ResultIterator(loadResultSet(queryString));
    }


    protected abstract T processResult(ResultSet rs) throws SQLException;

    private class ResultIterator implements CloseableTribbleIterator {

        private ResultSet rs;

        private boolean hasNext;
        private T next;

        public ResultIterator(ResultSet rs) {
            this.rs = rs;
            try {
                hasNext = rs.next();
            } catch (SQLException e) {
                log.error("Database error", e);
                throw new RuntimeException("Database error", e);
            }
        }

        @Override
        public boolean hasNext() {
            return hasNext;
        }

        @Override
        public T next() {
            try {
                next = processResult(rs);
                hasNext = rs.next();
                if (!hasNext) {
                    DBManager.closeResources(rs, rs.getStatement(), null);
                }

                return next;
            } catch (SQLException e) {
                log.error(e);
                throw new RuntimeException("Error processing SQL Results: " + e);
            }
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("Cannot remove");
        }

        @Override
        public Iterator<T> iterator() {
            return this;
        }

    }


}
