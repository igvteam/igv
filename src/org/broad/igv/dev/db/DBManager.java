package org.broad.igv.dev.db;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.LoginDialog;
import org.broad.igv.util.Utilities;
import org.w3c.dom.Attr;
import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.xml.sax.SAXException;

import javax.xml.parsers.ParserConfigurationException;
import java.awt.*;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.sql.*;
import java.util.*;
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
    static String username;
    static String password;

    private static Map<String, String> driverMap;
    static{
        driverMap = new HashMap<String, String>(2);
        driverMap.put("mysql", "com.mysql.jdbc.Driver");
        driverMap.put("sqlite", "org.sqlite.JDBC");
    }

    public static ConnectionWrapper getConnection(String url) {

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
        ConnectionWrapper conn = connect(url);
        if (conn != null) {
            connectionPool.put(url, conn);
            log.info("Connection pool size: " + connectionPool.size());
        }
        return conn;

    }

    static String createConnectionURL(String host, String db, String port, String subprotocol){
        String driver = driverMap.get(subprotocol);
        try {
            Class.forName(driver);
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }

        if(!host.startsWith("/")){
            host = "//" + host;
        }
        String url = "jdbc:" + subprotocol + ":" + host;
        if (port != null && !port.equals("")) {
            try{
                int iPort = Integer.parseInt(port);
                if(iPort >= 0){
                    url += ":" + iPort;
                }
            }catch(NumberFormatException e){
                log.error("Invalid port: " + port);
            }
        }
        if(db != null){
            url += "/" + db;
        }

        return url;
    }

    static Connection getConnection() {

        PreferenceManager preferenceManager = PreferenceManager.getInstance();
        String host = preferenceManager.get(PreferenceManager.DB_HOST);
        String db = preferenceManager.get(PreferenceManager.DB_NAME);
        String port = preferenceManager.get(PreferenceManager.DB_PORT);
        String subprotocol = "mysql";

        String url = createConnectionURL(host, db, port, subprotocol);
        return getConnection(url);
    }

    /**
     * Open connection using parameters specified in the given
     * profile.
     * @param profilePath
     * @return
     */
    public static Connection getStoredConnection(String profilePath){
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

            String url = createConnectionURL(host, path, port, subprotocol);
            return getConnection(url);

        } catch (Exception e) {
            e.printStackTrace();
            return null;
        } finally {
            try {
                if(profileStream != null) profileStream.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    private static String getNullSafe(NamedNodeMap attr, String key){
        Node node = attr.getNamedItem(key);
        return node != null ? node.getTextContent() : null;
    }

    private static ConnectionWrapper connect(String url) {
        try {
            return new ConnectionWrapper(DriverManager.getConnection(url, username, password));
        } catch (SQLException e) {
            int errorCode = e.getErrorCode();
            if (errorCode == 1044 || errorCode == 1045) {
                String host = PreferenceManager.getInstance().get(PreferenceManager.DB_HOST);

                Frame parent = Globals.isHeadless() ? null : IGV.getMainFrame();
                LoginDialog dlg = new LoginDialog(parent, false, host, false);
                dlg.setVisible(true);
                if (dlg.isCanceled()) {
                    throw new RuntimeException("Must login to access" + host);
                }
                username = dlg.getUsername();
                password = new String(dlg.getPassword());
                return connect(url);

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
