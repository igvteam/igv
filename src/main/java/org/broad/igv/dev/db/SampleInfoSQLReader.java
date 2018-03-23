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

    public SampleInfoSQLReader(DBProfile.DBTable table, String sampleColumn) {
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
