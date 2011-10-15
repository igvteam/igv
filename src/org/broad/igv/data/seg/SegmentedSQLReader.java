package org.broad.igv.data.seg;

import com.mysql.jdbc.Driver;
import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.LoginDialog;

import java.sql.*;

/**
 * Experimental class to explore using a SQL database as a data store
 *
 * @author Jim Robinson
 * @date 10/14/11
 */
public class SegmentedSQLReader {

    private static Logger log = Logger.getLogger(SegmentedSQLReader.class);

    public static void main(String[] args) {
        SegmentedAsciiDataSet ds = new SegmentedAsciiDataSet(null);
        (new SegmentedSQLReader()).load(ds);
    }


    public void load(SegmentedAsciiDataSet dataset) {


        System.out.println("MySQL Connect Example.");
        Connection conn = null;

        final PreferenceManager preferenceManager = PreferenceManager.getInstance();
        String host = preferenceManager.get(PreferenceManager.DB_HOST);
        String db = preferenceManager.get(PreferenceManager.DB_NAME);
        String port = preferenceManager.get(PreferenceManager.DB_PORT);

        String url = "jdbc:mysql://" + host;
        if (!port.equals("-1")) {
            url += ":" + port;
        }
        url += "/" + db;

        LoginDialog dlg = new LoginDialog(IGV.getMainFrame(), false, db, false);
        dlg.show();
        if (dlg.isCanceled()) {
            throw new RuntimeException("Must login to access" + db);
        }
        String userName = dlg.getUsername();
        String password = new String(dlg.getPassword());

        String driver = "com.mysql.jdbc.Driver";

        try {
            Class.forName(driver).newInstance();
            conn = DriverManager.getConnection(url, userName, password);
            System.out.println("Connected to the database");

            String query = "SELECT * FROM CNV";

            Statement st = conn.createStatement();
            ResultSet rs = st.executeQuery(query);
            while (rs.next()) {
                String sample = rs.getString("Sample");
                String chr = rs.getString("Chromosome");
                int start = Integer.parseInt(rs.getString("Start").replace(",", "")) - 1;
                int end = Integer.parseInt(rs.getString("Stop").replace(",", ""));
                float value = (float) Double.parseDouble(rs.getString("Probe Median"));
                String event = rs.getString("Event");
                String percentOverlap = rs.getString("% of CNV Overlap");
                String description = "<br>" +event + "<br>% CNV Overlap = " + percentOverlap;

                dataset.addSegment(sample, chr, start, end, value, description);
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
