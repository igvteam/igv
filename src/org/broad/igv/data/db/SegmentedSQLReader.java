package org.broad.igv.data.db;

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.data.seg.SegmentedAsciiDataSet;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.LoginDialog;
import org.broad.igv.util.ResourceLocator;

import java.sql.*;

/**
 * Experimental class to explore using a SQL database as a data store
 *
 * @author Jim Robinson
 * @date 10/14/11
 */
public class SegmentedSQLReader {

    private static Logger log = Logger.getLogger(SegmentedSQLReader.class);


    public SegmentedAsciiDataSet load(ResourceLocator locator, Genome genome) {

        System.out.println("MySQL Connect Example.");
        Connection conn = null;

        SegmentedAsciiDataSet dataset = new SegmentedAsciiDataSet(genome);

        final PreferenceManager preferenceManager = PreferenceManager.getInstance();
        String host = preferenceManager.get(PreferenceManager.DB_HOST);
        String db = preferenceManager.get(PreferenceManager.DB_NAME);
        String port = preferenceManager.get(PreferenceManager.DB_PORT);

        String url = locator.getServerURL();
        final String table = locator.getPath();

        LoginDialog dlg = new LoginDialog(IGV.getMainFrame(), false, db, false);
        dlg.show();
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
            while (rs.next()) {

                String sample = rs.getString("Sample");
                String chr = rs.getString("chr");
                int start = Integer.parseInt(rs.getString("start").replace(",", "")) - 1;
                int end = Integer.parseInt(rs.getString("end").replace(",", ""));
                float value = (float) Double.parseDouble(rs.getString("value"));
                String description =  rs.getString("description");
                dataset.addSegment(sample, chr, start, end, value, description);

            }
            dataset.sortLists();


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

        return dataset;
    }

}
