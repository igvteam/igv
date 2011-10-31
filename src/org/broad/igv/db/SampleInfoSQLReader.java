package org.broad.igv.db;

import org.apache.log4j.Logger;
import org.broad.igv.track.AttributeManager;

import org.broad.igv.util.ResourceLocator;

import java.sql.*;

/**
 * @author Jim Robinson
 * @date 10/15/11
 */
public class SampleInfoSQLReader {

    private static Logger log  = Logger.getLogger(SampleInfoSQLReader.class);

    public void load(ResourceLocator locator)  {


        ResultSet rs = null;
        Statement st = null;
        try {
            Connection conn = DBManager.getConnection();

            String query = locator.getDescription();
            st = conn.createStatement();
            rs = st.executeQuery(query);

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

       } catch (SQLException e) {
            log.error("Database error", e);
            throw new RuntimeException("Database error", e);
        } finally {
             if(rs != null) {
                 try {
                     rs.close();
                 } catch (SQLException e) {
                     e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                 }
             }
            if(st != null) {
                try {
                    st.close();
                } catch (SQLException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }
        }

    }
}
