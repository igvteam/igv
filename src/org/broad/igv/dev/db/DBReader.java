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
            colListing = StringUtils.join(colNames, ',');
        }
        return String.format("SELECT %s FROM %s", colListing, tableName);
    }

    public DBReader(DBTable table) {
        this(table.getDbLocator(), table.getTableName(), table.getBaseQuery());
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
            ResultSet rs = st.executeQuery(queryString);

            return rs;

        } catch (SQLException e) {
            log.error("Database error", e);
            throw new RuntimeException("Database error", e);
        }
    }

    public String getTableName() {
        return tableName;
    }

}
