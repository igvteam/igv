package org.broad.igv.db;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.LoginDialog;

import java.awt.*;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;

/**
 * Class for prototyping database connections.  Prototype only -- hardcoded for mysql,  connects to single database,
 * keeps single connection, etc.
 *
 * @author Jim Robinson
 * @date 10/31/11
 */
public class DBManager {

    private static Logger log = Logger.getLogger(DBManager.class);

    static Connection conn;
    static String username;
    static String password;




    public static Connection getConnection() {

        if (conn == null) {
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

            connect(url);
        }
        return conn;
    }

    private static void connect(String url) {
        try {
            conn = DriverManager.getConnection(url, username, password);
        } catch (SQLException e) {
            int errorCode = e.getErrorCode();
            if(errorCode == 1044 || errorCode == 1045) {
                String host = PreferenceManager.getInstance().get(PreferenceManager.DB_HOST);

                Frame parent = Globals.isHeadless() ? null : IGV.getMainFrame();
                LoginDialog dlg = new LoginDialog(parent, false, host, false);
                dlg.setVisible(true);
                if (dlg.isCanceled()) {
                    throw new RuntimeException("Must login to access" + host);
                }
                username = dlg.getUsername();
                password = new String(dlg.getPassword());
                connect(url);

            }
            else {
                MessageUtils.showMessage("<html>Error connecting to database: <br>" + e.getMessage());
            }

        }
    }


    public static void shutdown() {
        if (conn != null) {
            try {
                conn.close();
            } catch (SQLException e) {

            }
        }

    }


    /**
     * Test getting connection
     */

    public static void main(String[] args) {

        Globals.setHeadless(true);
        PreferenceManager preferenceManager = PreferenceManager.getInstance();
        preferenceManager.put(PreferenceManager.DB_HOST, "localhost");
        preferenceManager.put(PreferenceManager.DB_NAME, "nih");
        preferenceManager.put(PreferenceManager.DB_PORT, "-1");

        DBManager.getConnection();
        DBManager.shutdown();

    }

}
