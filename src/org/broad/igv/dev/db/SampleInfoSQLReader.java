package org.broad.igv.dev.db;

import org.apache.log4j.Logger;
import org.broad.igv.track.AttributeManager;
import org.broad.igv.util.ResourceLocator;

import java.sql.*;

/**
 * @author Jim Robinson
 * @date 10/15/11
 */
public class SampleInfoSQLReader extends DBReader {

    private static Logger log = Logger.getLogger(SampleInfoSQLReader.class);

    String sampleColumn = "SAMPLE_ID_ARRAY"; // TODO -- pass this in obviously

    public void load(ResourceLocator locator) {

        Connection conn = null;
        Statement st = null;
        ResultSet rs = null;

        try {
            conn = DBManager.getConnection();

            String query = locator.getDescription();
            st = conn.createStatement();
            rs = st.executeQuery(query);

            ResultSetMetaData metaData = rs.getMetaData();
            int nCol = metaData.getColumnCount();
            String[] columnNames = new String[nCol];
            for (int i = 0; i < nCol; i++) {
                columnNames[i] = metaData.getColumnName(i + 1);
            }

            while (rs.next()) {
                String sample = rs.getString(sampleColumn);
                for (String col : columnNames) {
                    String value = rs.getString(col);
                    AttributeManager.getInstance().addAttribute(sample, col, value);
                }
            }

        } catch (Exception e) {
            log.error("Error accessing database", e);
            throw new RuntimeException(e);
        } finally {
           closeResources(rs, st, conn);
        }

    }

}
