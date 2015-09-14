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
import org.broad.igv.util.ResourceLocator;

import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * Class for loading ALL of the data from a single table
 * in a single database.
 *
 * @author Jacob Silterra
 * @date 2012/05/30
 */
public abstract class WholeTableDBReader<T> extends DBReader {

    private static Logger log = Logger.getLogger(WholeTableDBReader.class);

    /**
     * See {@link DBReader#DBReader(DBProfile.DBTable)}
     */
    public WholeTableDBReader(DBProfile.DBTable table) {
        super(table);
    }

    /**
     * See {@link DBReader#DBReader(ResourceLocator, String, String)}
     */
    public WholeTableDBReader(ResourceLocator locator, String tableName, String baseQueryString) {
        super(locator, tableName, baseQueryString);
    }

    /**
     * Execute {@code #baseQueryString} and return the processed result
     * Connection is closed afterwards
     *
     * @return
     */
    public T load() {

        T obj = null;
        ResultSet rs = null;
        try {
            rs = super.executeQuery(baseQueryString);
            obj = processResultSet(rs);
        } catch (SQLException e) {
            log.error("Database error", e);
            throw new RuntimeException("Database error", e);
        } finally {
            DBManager.closeAll(rs);
        }

        return obj;
    }

    /**
     * Read off all relevant information from the ResultSet.
     * Closing the ResultSet is not necessary
     * @param rs
     * @return
     * @throws SQLException
     */
    protected abstract T processResultSet(ResultSet rs) throws SQLException;

}
