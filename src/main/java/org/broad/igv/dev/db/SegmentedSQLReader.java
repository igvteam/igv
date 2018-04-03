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
import org.broad.igv.data.seg.SegmentedAsciiDataSet;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.ResourceLocator;

import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * Reader for Segmented data, from a SQL database. Column names are hardcoded
 *
 * @author Jim Robinson
 * @date 10/14/11
 * @deprecated See {@link org.broad.igv.dev.SegmentedReader}
 */
public class SegmentedSQLReader extends WholeTableDBReader<SegmentedAsciiDataSet> {

    private static Logger log = Logger.getLogger(SegmentedSQLReader.class);
    private Genome genome;

    public SegmentedSQLReader(ResourceLocator locator, String tableName, Genome genome) {
        super(locator, tableName, null);
        this.genome = genome;
    }

    @Override
    protected SegmentedAsciiDataSet processResultSet(ResultSet rs) throws SQLException {
        SegmentedAsciiDataSet dataset = new SegmentedAsciiDataSet(genome);
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
        return dataset;
    }
}
