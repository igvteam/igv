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

import org.apache.log4j.Logger;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.StringUtils;

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Map;

/**
 * General class for reading from a SQL Database.
 * Each instance is attached to a specific database, table, and
 * configuration of columns.
 *
 * @author Jacob Silterra
 * @date 2012/05/30
 */
public class DBReader {

    private static Logger log = Logger.getLogger(DBReader.class);

    protected ResourceLocator locator;
    protected String baseQueryString;
    private String tableName;

    public static String buildBaseQueryString(String tableName, Map<Integer, String> columnLabelMap) {
        String colListing = "*";
        if (columnLabelMap != null) {
            String[] colNames = new String[columnLabelMap.size()];
            for (Map.Entry<Integer, String> entry : columnLabelMap.entrySet()) {
                colNames[entry.getKey()] = entry.getValue();
            }
            colListing = StringUtils.join(colNames, ",");
        }
        return String.format("SELECT %s FROM %s", colListing, tableName);
    }

    public DBReader(DBProfile.DBTable table) {
        this(table.getDbLocator(), table.getName(), table.getBaseQuery());
    }

    /**
     * @param locator
     * @param tableName
     * @param baseQueryString Can be null, in which case all columns from {@code tableName} are selected
     */
    public DBReader(ResourceLocator locator, String tableName, String baseQueryString) {
        this.locator = locator;
        assert tableName != null;
        this.tableName = tableName;
        this.baseQueryString = baseQueryString;
        if (baseQueryString == null) {
            this.baseQueryString = buildBaseQueryString(this.tableName, null);
        }
    }

    protected ResultSet executeQuery(String queryString) {

        try {
            Connection conn = DBManager.getConnection(locator);
            Statement st = conn.createStatement();
            return st.executeQuery(queryString);
        } catch (SQLException e) {
            log.error("Database error", e);
            throw new RuntimeException("Database error", e);
        }
    }

    public String getTableName() {
        return tableName;
    }

}
