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
