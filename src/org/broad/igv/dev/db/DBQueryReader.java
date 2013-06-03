/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.dev.db;

import org.apache.log4j.Logger;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Class for reading only portions of a table (queries) repeatedly.
 * The connection is NOT closed between queries.
 *
 * @author Jacob Silterra
 * @date 2012/05/30
 */
public abstract class DBQueryReader<T> extends DBReader {

    private static Logger log = Logger.getLogger(DBQueryReader.class);
    protected DBProfile.DBTable table;

    public DBQueryReader(DBProfile.DBTable table) {
        super(table);
        this.table = table;
    }

    protected Iterator loadIterator(PreparedStatement st) {
        try {
            ResultSet rs = st.executeQuery();
            return loadIterator(rs);
        } catch (SQLException e) {
            log.error("Database error", e);
            throw new RuntimeException("Database error", e);
        }
    }

    /**
     * Load all results from the provided {@code ResultSet}
     * so we can close it quickly
     * @param rs
     * @return
     */
    protected Iterator loadIterator(ResultSet rs) {
        List<T> results = new ArrayList<T>();
        try {
            while (rs.next()) {
                results.add(processResult(rs));
            }
        } catch (SQLException e) {
            log.error("Database error", e);
            throw new RuntimeException("Database error", e);
        } finally {
            DBManager.closeAll(rs);
        }
        return results.iterator();
    }

    /**
     * Read a single line from the {@code ResultSet} and return the object
     * based on it.
     * @param rs
     * @return
     * @throws SQLException
     */
    protected abstract T processResult(ResultSet rs) throws SQLException;

}
