package org.broad.igv.dev.db;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.LoginDialog;

import java.awt.*;
import java.sql.*;
import java.util.*;

/**
 * Class for prototyping database connections.  Prototype only -- hardcoded for mysql,  connects to single database,
 * keeps single connection, etc.
 *
 * @author Jim Robinson
 * @date 10/31/11
 */
public class DBManager {

    private static Logger log = Logger.getLogger(DBManager.class);

    static Collection<ConnectionWrapper> connectionPool =
            Collections.synchronizedCollection(new ArrayList<ConnectionWrapper>());
    static String username;
    static String password;

    public static Connection getConnection() {

        Iterator<ConnectionWrapper> poolIter = connectionPool.iterator();
        while (poolIter.hasNext()) {
            ConnectionWrapper conn = poolIter.next();
            try {
                if (conn == null || conn.isReallyClosed()) {
                    poolIter.remove();
                } else if (!conn.isClosed()) {
                    return conn;
                }
            } catch (SQLException e) {
                log.error("Bad connection", e);
                poolIter.remove();
            }
        }

        // No valid connections
        ConnectionWrapper conn = createConnection();
        if (conn != null) {
            connectionPool.add(conn);
            log.info("Connection pool size: " + connectionPool.size());
        }
        return conn;

    }


    private static ConnectionWrapper createConnection() {
        String driver = "com.mysql.jdbc.Driver";
        try {
            Class.forName(driver).newInstance();
        } catch (Exception e) {
            e.printStackTrace();
        }

        final PreferenceManager preferenceManager = PreferenceManager.getInstance();
        String host = preferenceManager.get(PreferenceManager.DB_HOST);
        String db = preferenceManager.get(PreferenceManager.DB_NAME);
        String port = preferenceManager.get(PreferenceManager.DB_PORT);

        String url = "jdbc:mysql://" + host;
        if (!port.equals("-1")) {
            url += ":" + port;
        }
        url += "/" + db;

        return connect(url);
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
        for (ConnectionWrapper conn : connectionPool) {
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
