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
     * See {@link DBReader#DBReader(DBTable)}
     */
    public WholeTableDBReader(DBTable table) {
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
            DBManager.closeResources(rs, null, DBManager.getConnection(locator));
        }

        return obj;
    }


    protected abstract T processResultSet(ResultSet rs) throws SQLException;

}
