package org.broad.igv.data.db;

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.track.AttributeManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.LoginDialog;
import org.broad.igv.util.ResourceLocator;

import java.sql.*;

/**
 * @author Jim Robinson
 * @date 10/15/11
 */
public class SampleInfoSQLReader {

    private static Logger log  = Logger.getLogger(SampleInfoSQLReader.class);

    public void load(ResourceLocator locator) {

        System.out.println("MySQL Connect Example.");
        Connection conn = null;


        final PreferenceManager preferenceManager = PreferenceManager.getInstance();
        String db = preferenceManager.get(PreferenceManager.DB_NAME);

        String url = locator.getServerURL();
        final String table = locator.getPath();

        LoginDialog dlg = new LoginDialog(IGV.getMainFrame(), false, db, false);
        dlg.setVisible(true);
        if (dlg.isCanceled()) {
            throw new RuntimeException("Must login to access" + db);
        }
        String userName = dlg.getUsername();
        String password = new String(dlg.getPassword());

        // TODO -- look based on url, or let user provide this
        String driver = "com.mysql.jdbc.Driver";

        try {
            Class.forName(driver).newInstance();
            conn = DriverManager.getConnection(url, userName, password);
            System.out.println("Connected to the database");

            String query = locator.getDescription();

            Statement st = conn.createStatement();
            ResultSet rs = st.executeQuery(query);

            ResultSetMetaData metaData = rs.getMetaData();
            int nCol = metaData.getColumnCount();
            String [] columnNames = new String[nCol];
            for(int i=0; i<nCol; i++) {
                columnNames[i] = metaData.getColumnName(i+1);
            }

            while (rs.next()) {

                String sample = rs.getString("Sample");
                for(String col : columnNames) {
                    String value = rs.getString(col);
                    AttributeManager.getInstance().addAttribute(sample, col, value);
                }


            }


            System.out.println("Disconnected from database");
        } catch (Exception e) {
            log.error("Error accessing database", e);
            throw new RuntimeException(e);
        } finally {
            if (conn != null) try {
                conn.close();
            } catch (SQLException e) {
                log.error("Error closing sql connection", e);
            }
        }

    }
}
