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
import org.broad.igv.Globals;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.LoginDialog;
import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.io.File;
import java.sql.*;
import java.util.*;

/**
 * Class for database connections.  Manages connections based on
 * ResourceLocators, connects, loads drivers, and can perform limited data conversion (from ResultSet -> String, String[])
 *
 * @author Jim Robinson
 * @date 10/31/11
 */
public class DBManager {

    private static Logger log = Logger.getLogger(DBManager.class);

    static Map<String, Connection> connectionPool =
            Collections.synchronizedMap(new HashMap<String, Connection>());

    private static Map<String, String> driverMap;

    static {
        driverMap = new HashMap<String, String>(2);
        driverMap.put("mysql", "com.mysql.jdbc.Driver");
        driverMap.put("sqlite", "org.sqlite.JDBC");
        driverMap.put("oracle:thin", "oracle.jdbc.driver.OracleDriver");
        driverMap.put("oracle", "oracle.jdbc.driver.OracleDriver");
    }

    public static Connection getConnection(ResourceLocator locator) {
        String url = locator.getPath();
        if (connectionPool.containsKey(url)) {
            Connection conn = connectionPool.get(url);
            try {
                if (conn == null || conn.isClosed()) {
                    connectionPool.remove(url);
                } else {
                    return conn;
                }
            } catch (SQLException e) {
                log.error("Bad connection", e);
                connectionPool.remove(url);
            }
        }

        // No valid connections
        Connection conn = connect(locator);
        if (conn != null) {
            connectionPool.put(url, conn);
            log.info("Connection pool size: " + connectionPool.size());
        }
        return conn;

    }

    public static void closeConnection(ResourceLocator locator) {
        String url = locator.getPath();
        if (connectionPool.containsKey(url)) {
            Connection conn = connectionPool.get(url);
            try {
                if (conn != null && !conn.isClosed()) {
                    conn.close();
                    connectionPool.remove(url);
                }
            } catch (SQLException e) {
                log.error(e.getMessage(), e);
            }
        }
    }

    private static String getSubprotocol(String url) {
        String[] tokens = url.split(":");
        return tokens[1];
    }

    private static Connection connect(ResourceLocator locator) {
        createDriver(getSubprotocol(locator.getPath()));
        try {
            return DriverManager.getConnection(locator.getPath(),
                    locator.getUsername(), locator.getPassword());
        } catch (SQLException e) {
            int errorCode = e.getErrorCode();
            if (errorCode == 1044 || errorCode == 1045) {
                String resource = locator.getPath();

                Frame parent = Globals.isHeadless() ? null : IGV.getMainFrame();
                LoginDialog dlg = new LoginDialog(parent, false, resource, false);
                dlg.setVisible(true);
                if (dlg.isCanceled()) {
                    throw new RuntimeException("Must login to access" + resource);
                }
                locator.setUsername(dlg.getUsername());
                locator.setPassword(new String(dlg.getPassword()));
                return connect(locator);

            } else {
                MessageUtils.showErrorMessage("<html>Error connecting to database: <br>" + e.getMessage(), e);
                return null;
            }

        }
    }

    public static void shutdown() {
        for (Connection conn : connectionPool.values()) {
            if (conn != null) {
                try {
                    conn.close();
                } catch (SQLException e) {
                    log.error("Error shutting down " + conn.toString(), e);
                }
            }
        }
        connectionPool.clear();
    }

    public static java.lang.Class<?> createDriver(String subprotocol) {
        String driver = driverMap.get(subprotocol);
        try {
            return Class.forName(driver);
        } catch (ClassNotFoundException e) {
            log.error("Unable to create driver for " + subprotocol, e);
            throw new IllegalArgumentException(e);
        }
    }


    /**
     * Creates a connection URL, and loads the driver class necessary for the protocol
     *
     * @param subprotocol
     * @param host
     * @param db
     * @param port
     * @return
     */
    public static String createConnectionURL(String subprotocol, String host, String db, String port) {
        createDriver(subprotocol);

        //If the host is a local file or oracle, don't want the leading "//"
        boolean isOracle = subprotocol.toLowerCase().contains("oracle");
        if (!isOracle && !(new File(host)).exists()) {
            host = "//" + host;
        }
        if(isOracle){
            host = "@" + subprotocol;
        }
        String url = "jdbc:" + subprotocol + ":" + host;
        if (port != null && !port.equals("")) {
            try {
                int iPort = Integer.parseInt(port);
                if (iPort >= 0) {
                    url += ":" + iPort;
                }
            } catch (NumberFormatException e) {
                log.error("Invalid port: " + port, e);
            }
        }
        if (db != null) {
            url += "/" + db;
        }

        return url;
    }

    /**
     * Close all resources associated with the ResultSet,
     * including the statement and connection
     * @param rs
     */
    public static void closeAll(ResultSet rs){
        Statement st = null;
        Connection conn = null;

        if (rs != null) {
            try {
                st = rs.getStatement();
                conn = st.getConnection();
            } catch (SQLException e) {
                log.error("Error getting statement and connection from result set", e);
            }
        }

        DBManager.closeResources(rs, st, conn);
    }

    /**
     * Close the specified resources
     *
     * @param rs
     * @param st
     * @param conn
     */
    static void closeResources(ResultSet rs, Statement st, Connection conn) {
        if (rs != null) {
            try {
                rs.close();
            } catch (SQLException e) {
                log.error("Error closing resultset", e);
            }
        }
        if (st != null) {
            try {
                st.close();
            } catch (SQLException e) {
                log.error("Error closing statement", e);
            }
        }
        if (conn != null) {
            try {
                conn.close();
                pruneConnectionPool();
            } catch (SQLException e) {
                log.error("Error closing sql connection", e);
            }
        }

    }

    /**
     * Remove all closed connection objects from the connection pool.
     * @throws SQLException
     */
    private static void pruneConnectionPool() throws SQLException{
        Set<String> toRemove = new HashSet<String>();
        for(Map.Entry<String, Connection> entry: connectionPool.entrySet()){
            if(entry.getValue().isClosed()){
                toRemove.add(entry.getKey());
            }
        }
        for(String url: toRemove){
            connectionPool.remove(url);
        }
    }

    /**
     * Convert the ResultSet into a string array, re-arranging columns according to
     * {@code columnIndexMap}, which is a map from array indexes -> sql column indexes
     *
     * @param rs
     * @param columnIndexMap
     * @return
     * @throws SQLException
     */
    public static String[] lineToArray(ResultSet rs, Map<Integer, String> columnIndexMap) throws SQLException {
        String[] colNames = DBProfile.DBTable.columnMapToArray(columnIndexMap);
        String[] tokens = new String[colNames.length];

        for (int cc = 0; cc < colNames.length; cc++) {
            String sqlCol = colNames[cc];
            if (sqlCol != null) {
                tokens[cc] = getStringFromResultSet(rs, sqlCol);
            }
        }
        return tokens;
    }

    /**
     * Convert a the current line to an array of strings
     *
     * @param rs
     * @param startColIndex 1-based start column index (lower columns are ignored)
     * @param endColIndex   1-based, inclusive end column index (columns afterwards are ignored)
     * @param labelsOnly    Whether to only include the column labels, rather than the data contained in the ResultSet
     * @return
     * @throws SQLException
     */
    public static String[] lineToArray(ResultSet rs, int startColIndex, int endColIndex, boolean labelsOnly) throws SQLException {
        ResultSetMetaData md = rs.getMetaData();
        int colCount = Math.min(md.getColumnCount(), endColIndex) - startColIndex + 1;
        String[] tokens = new String[colCount];
        for (int cc = 0; cc < colCount; cc++) {
            int sqlCol = cc + startColIndex;
            if (labelsOnly) {
                tokens[cc] = md.getColumnLabel(sqlCol);
            } else {
                tokens[cc] = getStringFromResultSet(rs, sqlCol);
            }
        }
        return tokens;
    }

    private static String getStringFromResultSet(ResultSet rs, String columnLabel) throws SQLException {
        return getStringFromResultSet(rs, rs.findColumn(columnLabel));
    }


    /**
     * Get the value at column {@code sqlCol} in the current row as a string.
     * <p/>
     * Have to parse blobs specially, otherwise we get the pointer as a string
     *
     * @param rs
     * @param sqlCol 1-indexed column number
     * @return
     * @throws SQLException
     */
    private static String getStringFromResultSet(ResultSet rs, int sqlCol) throws SQLException {
        String s;
        int type = rs.getMetaData().getColumnType(sqlCol);

        if (blobTypes.contains(type)) {
            Blob b = rs.getBlob(sqlCol);
            s = new String(b.getBytes(1l, (int) b.length()));
        } else {
            s = rs.getString(sqlCol);
        }
        return s;
    }

    private static final Set<Integer> blobTypes;

    static {
        int[] blobtypes = {Types.BINARY, Types.BLOB, Types.VARBINARY, Types.LONGVARBINARY};
        blobTypes = new HashSet<Integer>(blobtypes.length);
        for (int bt : blobtypes) {
            blobTypes.add(bt);
        }
    }

}
