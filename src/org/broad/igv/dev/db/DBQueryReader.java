/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
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
