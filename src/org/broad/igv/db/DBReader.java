package org.broad.igv.db;

import org.apache.log4j.Logger;

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

/**
 * @author Jim Robinson
 * @date 1/20/12
 */
public abstract class DBReader {

    private static Logger log = Logger.getLogger(DBReader.class);


    public static void closeResources(ResultSet rs, Statement st, Connection conn) {
        if (rs != null) {
            try {
                rs.close();
            } catch (SQLException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }
        if (st != null) {
            try {
                st.close();
            } catch (SQLException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }
        if (conn != null) try {
            conn.close();
        } catch (SQLException e) {
            log.error("Error closing sql connection", e);
        }

    }
}
