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
import org.broad.igv.Globals;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.LoginDialog;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.Utilities;
import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;

import java.awt.*;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.sql.*;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Properties;
import java.util.concurrent.Executor;

/**
 * Class for prototyping database connections.  Prototype only -- hardcoded for mysql,  connects to single database,
 * keeps single connection, etc.
 *
 * @author Jim Robinson
 * @date 10/31/11
 */
public class DBManager {

    private static Logger log = Logger.getLogger(DBManager.class);

    static Map<String, ConnectionWrapper> connectionPool =
            Collections.synchronizedMap(new HashMap<String, ConnectionWrapper>());

    private static Map<String, String> driverMap;

    static {
        driverMap = new HashMap<String, String>(2);
        driverMap.put("mysql", "com.mysql.jdbc.Driver");
        driverMap.put("sqlite", "org.sqlite.JDBC");
    }

    public static ConnectionWrapper getConnection(ResourceLocator locator) {
        String url = locator.getPath();
        if (connectionPool.containsKey(url)) {
            ConnectionWrapper conn = connectionPool.get(url);
            try {
                if (conn == null || conn.isReallyClosed()) {
                    connectionPool.remove(url);
                } else if (!conn.isClosed()) {
                    return conn;
                }
            } catch (SQLException e) {
                log.error("Bad connection", e);
                connectionPool.remove(url);
            }
        }


        // No valid connections
        ConnectionWrapper conn = connect(locator);
        if (conn != null) {
            connectionPool.put(url, conn);
            log.info("Connection pool size: " + connectionPool.size());
        }
        return conn;

    }

    private static String getNullSafe(NamedNodeMap attr, String key) {
        Node node = attr.getNamedItem(key);
        return node != null ? node.getTextContent() : null;
    }

    private static ConnectionWrapper connect(ResourceLocator locator) {
        try {
            return new ConnectionWrapper(DriverManager.getConnection(locator.getPath(),
                    locator.getUsername(), locator.getPassword()));
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
                MessageUtils.showMessage("<html>Error connecting to database: <br>" + e.getMessage());
                return null;
            }

        }
    }

    public static void shutdown() {
        for (ConnectionWrapper conn : connectionPool.values()) {
            if (conn != null) {
                try {
                    conn.reallyClose();
                } catch (SQLException e) {

                }
            }
        }
        connectionPool.clear();
    }


    static String createConnectionURL(String subprotocol, String host, String db, String port) {
        String driver = driverMap.get(subprotocol);
        try {
            Class.forName(driver);
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }

        //If the host is a local file, don't want the leading "//"
        if (!(new File(host)).exists()) {
            host = "//" + host;
        }
        String url = "jdbc:" + subprotocol + ":" + host;
        if (port != null && !port.equals("")) {
            try {
                int iPort = Integer.parseInt(port);
                if (iPort >= 0) {
                    url += ":" + iPort;
                }
            } catch (NumberFormatException e) {
                log.error("Invalid port: " + port);
            }
        }
        if (db != null) {
            url += "/" + db;
        }

        return url;
    }

    /**
     * Open connection using parameters specified in the given
     * profile.
     *
     * @param profilePath
     * @return
     */
    public static ResourceLocator getStoredConnection(String profilePath) {
        InputStream profileStream = null;
        try {
            profileStream = new FileInputStream(profilePath);
            Document document = Utilities.createDOMDocumentFromXmlStream(profileStream);
            Node db = document.getElementsByTagName("database").item(0);
            NamedNodeMap attr = db.getAttributes();
            String host = attr.getNamedItem("host").getTextContent();
            String path = attr.getNamedItem("path").getTextContent();
            String subprotocol = attr.getNamedItem("subprotocol").getTextContent();

            String port = getNullSafe(attr, "port");
            String username = getNullSafe(attr, "username");
            String password = getNullSafe(attr, "password");

            ResourceLocator locator = new ResourceLocator(createConnectionURL(subprotocol, host, path, port));
            locator.setUsername(username);
            locator.setPassword(password);

            return locator;

        } catch (Exception e) {
            e.printStackTrace();
            return null;
        } finally {
            try {
                if (profileStream != null) profileStream.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }


    static class ConnectionWrapper implements Connection {

        Connection conn;
        boolean closed;

        ConnectionWrapper(Connection conn) {
            this.conn = conn;
            closed = false;
        }

        public void close() throws SQLException {
            closed = true;
        }

        public boolean isClosed() throws SQLException {
            return closed;
        }

        public void reallyClose() throws SQLException {
            closed = true;
            conn.close();
        }

        public boolean isReallyClosed() throws SQLException {
            return conn.isClosed();
        }


        public void clearWarnings() throws SQLException {
            conn.clearWarnings();
        }

        public void commit() throws SQLException {
            conn.commit();
        }

        public Array createArrayOf(String s, Object[] objects) throws SQLException {
            return conn.createArrayOf(s, objects);
        }

        public Blob createBlob() throws SQLException {
            return conn.createBlob();
        }

        public Clob createClob() throws SQLException {
            return conn.createClob();
        }

        public NClob createNClob() throws SQLException {
            return conn.createNClob();
        }

        public SQLXML createSQLXML() throws SQLException {
            return conn.createSQLXML();
        }

        public Statement createStatement() throws SQLException {
            return conn.createStatement();
        }

        public Statement createStatement(int i, int i1) throws SQLException {
            return conn.createStatement(i, i1);
        }

        public Statement createStatement(int i, int i1, int i2) throws SQLException {
            return conn.createStatement(i, i1, i2);
        }

        public Struct createStruct(String s, Object[] objects) throws SQLException {
            return conn.createStruct(s, objects);
        }

        public void setSchema(String schema) throws SQLException {
            throw new UnsupportedOperationException("Operation not supported for backwards compatibility to Java 6");
        }

        public String getSchema() throws SQLException {
            return null; //TODO
        }

        public void abort(Executor executor) throws SQLException {
            throw new UnsupportedOperationException("Operation not supported for backwards compatibility to Java 6");
        }

        public void setNetworkTimeout(Executor executor, int milliseconds) throws SQLException {
            throw new UnsupportedOperationException("Operation not supported for backwards compatibility to Java 6");
        }

        public int getNetworkTimeout() {
            return -1;
        }

        public boolean getAutoCommit() throws SQLException {
            return conn.getAutoCommit();
        }

        public String getCatalog() throws SQLException {
            return conn.getCatalog();
        }

        public Properties getClientInfo() throws SQLException {
            return conn.getClientInfo();
        }

        public String getClientInfo(String s) throws SQLException {
            return conn.getClientInfo(s);
        }

        public int getHoldability() throws SQLException {
            return conn.getHoldability();
        }

        public DatabaseMetaData getMetaData() throws SQLException {
            return conn.getMetaData();
        }

        public int getTransactionIsolation() throws SQLException {
            return conn.getTransactionIsolation();
        }

        public Map<String, Class<?>> getTypeMap() throws SQLException {
            return conn.getTypeMap();
        }

        public SQLWarning getWarnings() throws SQLException {
            return conn.getWarnings();
        }

        public boolean isReadOnly() throws SQLException {
            return conn.isReadOnly();
        }

        public boolean isValid(int i) throws SQLException {
            return conn.isValid(i);
        }

        public String nativeSQL(String s) throws SQLException {
            return conn.nativeSQL(s);
        }

        public CallableStatement prepareCall(String s) throws SQLException {
            return conn.prepareCall(s);
        }

        public CallableStatement prepareCall(String s, int i, int i1) throws SQLException {
            return conn.prepareCall(s, i, i1);
        }

        public CallableStatement prepareCall(String s, int i, int i1, int i2) throws SQLException {
            return conn.prepareCall(s, i, i1, i2);
        }

        public PreparedStatement prepareStatement(String s) throws SQLException {
            return conn.prepareStatement(s);
        }

        public PreparedStatement prepareStatement(String s, int i) throws SQLException {
            return conn.prepareStatement(s, i);
        }

        public PreparedStatement prepareStatement(String s, int i, int i1) throws SQLException {
            return conn.prepareStatement(s, i, i1);
        }

        public PreparedStatement prepareStatement(String s, int i, int i1, int i2) throws SQLException {
            return conn.prepareStatement(s, i, i1, i2);
        }

        public PreparedStatement prepareStatement(String s, int[] ints) throws SQLException {
            return conn.prepareStatement(s, ints);
        }

        public PreparedStatement prepareStatement(String s, String[] strings) throws SQLException {
            return conn.prepareStatement(s, strings);
        }

        public void releaseSavepoint(Savepoint savepoint) throws SQLException {
            conn.releaseSavepoint(savepoint);
        }

        public void rollback() throws SQLException {
            conn.rollback();
        }

        public void rollback(Savepoint savepoint) throws SQLException {
            conn.rollback(savepoint);
        }

        public void setAutoCommit(boolean b) throws SQLException {
            conn.setAutoCommit(b);
        }

        public void setCatalog(String s) throws SQLException {
            conn.setCatalog(s);
        }

        public void setClientInfo(Properties properties) throws SQLClientInfoException {
            conn.setClientInfo(properties);
        }

        public void setClientInfo(String s, String s1) throws SQLClientInfoException {
            conn.setClientInfo(s, s1);
        }

        public void setHoldability(int i) throws SQLException {
            conn.setHoldability(i);
        }

        public void setReadOnly(boolean b) throws SQLException {
            conn.setReadOnly(b);
        }

        public Savepoint setSavepoint() throws SQLException {
            return conn.setSavepoint();
        }

        public Savepoint setSavepoint(String s) throws SQLException {
            return conn.setSavepoint(s);
        }

        public void setTransactionIsolation(int i) throws SQLException {
            conn.setTransactionIsolation(i);
        }

        public void setTypeMap(Map<String, Class<?>> stringClassMap) throws SQLException {
            conn.setTypeMap(stringClassMap);
        }

        public boolean isWrapperFor(Class<?> aClass) throws SQLException {
            return conn.isWrapperFor(aClass);
        }

        public <T> T unwrap(Class<T> tClass) throws SQLException {
            return conn.unwrap(tClass);
        }

    }


}
