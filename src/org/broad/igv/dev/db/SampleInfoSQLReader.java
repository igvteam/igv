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
import org.broad.igv.track.AttributeManager;
import org.broad.igv.util.ResourceLocator;

import java.sql.ResultSet;
import java.sql.ResultSetMetaData;
import java.sql.SQLException;

/**
 * @author Jim Robinson
 * @date 10/15/11
 */
public class SampleInfoSQLReader extends WholeTableDBReader<Void> {

    private static Logger log = Logger.getLogger(SampleInfoSQLReader.class);

    String sampleColumn; // = "SAMPLE_ID_ARRAY";

    public SampleInfoSQLReader(DBTable table, String sampleColumn) {
        super(table);
        this.sampleColumn = sampleColumn;
    }

    public SampleInfoSQLReader(ResourceLocator locator, String tableName, String sampleColumn) {
        super(locator, tableName, null);
        this.sampleColumn = sampleColumn;
    }

    @Override
    protected Void processResultSet(ResultSet rs) throws SQLException {
        ResultSetMetaData metaData = rs.getMetaData();
        int nCol = metaData.getColumnCount();
        String[] columnNames = new String[nCol];
        for (int i = 0; i < nCol; i++) {
            columnNames[i] = metaData.getColumnName(i + 1);
        }

        while (rs.next()) {
            String sample = rs.getString(sampleColumn);
            for (String col : columnNames) {
                if (col.equals(sampleColumn)) {
                    continue;
                }
                String value = rs.getString(col);
                AttributeManager.getInstance().addAttribute(sample, col, value);
            }
        }

        return null;
    }
}
