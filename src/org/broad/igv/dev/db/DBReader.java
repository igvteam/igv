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

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.CloseableTribbleIterator;

import java.sql.*;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

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
    protected String tableName;
    protected String baseQueryString = "SELECT * FROM ";
    protected ColumnMap columnMap;

    public DBReader(ResourceLocator locator, String tableName, ColumnMap columnMap) {
        this.locator = locator;
        assert tableName != null;
        this.tableName = tableName;
        this.columnMap = columnMap;
        String colListing = "*";
        if (this.columnMap != null) {
            String[] colNames = new String[columnMap.columnLabelMap.size()];
            for (Map.Entry<Integer, String> entry : columnMap.columnLabelMap.entrySet()) {
                colNames[entry.getKey()] = entry.getValue();
            }
            colListing = StringUtils.join(colNames, ',');
        }
        this.baseQueryString = String.format("SELECT %s FROM %s", colListing, this.tableName);
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
        return new ResultIterator(loadResultSet(queryString));
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

    public String getTableName() {
        return tableName;
    }

    protected final int getDBColumn(int fileColNum) {
        if (this.columnMap != null) {
            return this.columnMap.getDBColumn(fileColNum);
        }
        return fileColNum;
    }

    public static class ColumnMap {
        private Map<Integer, Integer> columnIndexMap = new HashMap<Integer, Integer>();
        private Map<Integer, String> columnLabelMap = new HashMap<Integer, String>();
        int minFileColNum = Integer.MAX_VALUE;
        int maxFileColNum = -1;

        private void put(int fileColNum, int dbColNum) {
            columnIndexMap.put(fileColNum, dbColNum);
            if (dbColNum < minFileColNum) minFileColNum = fileColNum;
            if (dbColNum > maxFileColNum) maxFileColNum = fileColNum;
        }

        void put(int fileColNum, String dbLabel) {
            if (dbLabel != null) {
                columnLabelMap.put(fileColNum, dbLabel);
            }
        }

        public void labelsToIndexes(ResultSetMetaData metaData) throws SQLException {
            for (Map.Entry<Integer, String> labelEntry : columnLabelMap.entrySet()) {
                String label = labelEntry.getValue();
                int index = findColumnByLabel(metaData, label);
                if (index < 0) {
                    throw new SQLException("Column " + label + " not found");
                }
                put(labelEntry.getKey(), index);
            }
        }

        private int findColumnByLabel(ResultSetMetaData metaData, String label) throws SQLException {
            for (int cc = 1; cc <= metaData.getColumnCount(); cc++) {
                if (metaData.getColumnLabel(cc).equals(label)) {
                    return cc;
                }
            }
            return -1;
        }

        public int getDBColumn(int fileColNum) {
            return columnIndexMap.get(fileColNum);
        }

        public int getMaxFileColNum() {
            return maxFileColNum;
        }

        public int getMinFileColNum() {
            return minFileColNum;
        }

    }


}
