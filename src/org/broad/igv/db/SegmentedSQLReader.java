package org.broad.igv.db;

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
public class SegmentedSQLReader extends DBReader {

    private static Logger log = Logger.getLogger(SegmentedSQLReader.class);


    public SegmentedAsciiDataSet load(ResourceLocator locator, Genome genome) {

        SegmentedAsciiDataSet dataset = new SegmentedAsciiDataSet(genome);

        ResultSet rs = null;
        Statement st = null;
        Connection conn = null;
        try {
            conn = DBManager.getConnection();
            String query = locator.getDescription();
            st = conn.createStatement();
            rs = st.executeQuery(query);
            while (rs.next()) {

                String sample = rs.getString("Sample");
                String chr = rs.getString("chr");
                int start = rs.getInt("start");
                int end = rs.getInt("end");
                float value = rs.getFloat("value");
                String description = rs.getString("description");
                dataset.addSegment(sample, chr, start, end, value, description);

            }
            dataset.sortLists();


            System.out.println("Disconnected from database");
        } catch (SQLException e) {
            log.error("Database error", e);
            throw new RuntimeException("Database error", e);
        } finally {
            closeResources(rs, st, conn);

        }

        return dataset;
    }

}
