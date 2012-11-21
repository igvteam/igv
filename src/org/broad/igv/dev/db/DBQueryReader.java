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
import org.broad.tribble.CloseableTribbleIterator;

import java.sql.*;
import java.util.Iterator;

/**
 * Class for reading only portions of a table (queries) repeatedly.
 * The connection is NOT closed between queries.
 *
 * @author Jacob Silterra
 * @date 2012/05/30
 */
public abstract class DBQueryReader<T> extends DBReader {

    private static Logger log = Logger.getLogger(DBQueryReader.class);
    protected DBTable table;

    public DBQueryReader(DBTable table) {
        super(table);
        this.table = table;
    }

    protected ResultSet executeQuery(String queryString) {

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

    protected CloseableTribbleIterator loadIterator(PreparedStatement st) {
        try {
            ResultSet rs = st.executeQuery();

            return new ResultIterator(rs, false);

        } catch (SQLException e) {
            log.error("Database error", e);
            throw new RuntimeException("Database error", e);
        }
    }

    protected CloseableTribbleIterator loadIterator(String queryString) {
        return new ResultIterator(executeQuery(queryString));
    }

    protected abstract T processResult(ResultSet rs) throws SQLException;

    private class ResultIterator implements CloseableTribbleIterator {

        private ResultSet rs;
        /**
         * Whether the statement will be closed when iteration
         * is complete. The ResultSet is always closed
         */
        private boolean closeStatement;

        private boolean hasNext;
        private T next;

        public ResultIterator(ResultSet rs) {
            this(rs, true);
        }

        public ResultIterator(ResultSet rs, boolean closeStatement) {
            this.rs = rs;
            this.closeStatement = closeStatement;
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
                    DBManager.closeResources(rs, closeStatement ? rs.getStatement() : null, null);
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

        @Override
        public void close() {
            try {
                DBManager.closeResources(rs, rs.getStatement(), null);
            } catch (SQLException e) {
                log.error(e);
                throw new RuntimeException(e);
            }
        }
    }


}
