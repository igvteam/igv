package org.broad.igv.db;

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

    String sampleColumn = "SAMPLE_ID_ARRAY"; // TODO -- pass this in obviously

    public void load(ResourceLocator locator) {

        System.out.println("MySQL Connect Example.");
        Connection conn = null;

        try {
            conn = DBManager.getConnection();

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

                String sample = rs.getString(sampleColumn);
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
