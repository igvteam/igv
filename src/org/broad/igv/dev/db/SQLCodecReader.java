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

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

/**
 * Class for reading data from SQL database, where the lines
 * are of a format we can read with existing codecs
 *
 * @author Jacob Silterra
 * @date 29 May 2012
 */
public class SQLCodecReader extends DBReader<Iterable<Feature>> {

    private static Logger log = Logger.getLogger(SQLCodecReader.class);

    protected FeatureCodec codec;

    public SQLCodecReader(FeatureCodec codec) {
        this.codec = codec;
    }

    @Override
    protected Iterable<Feature> processResultSet(ResultSet rs) throws SQLException {

        List<Feature> featureList;
        featureList = new ArrayList<Feature>(rs.getFetchSize());
        while (rs.next()) {
            Feature feat = this.parseLine(rs);
            featureList.add(feat);
        }

        return featureList;
    }

    /**
     * Turn a line from a result set into a feature
     * TODO Make abstract
     *
     * @param rs
     * @return
     * @throws SQLException
     */
    protected Feature parseLine(ResultSet rs) throws SQLException {
        String[] tokens = lineToArray(rs);
        //TODO GET RID OF THIS, IT'S BAD AND I FEEL BAD FOR WRITING IT -JS
        String line = StringUtils.join(tokens, "\t");
        return codec.decode(line);

    }

    /**
     * Convert a the current line to an array of strings
     *
     * @param rs
     * @return
     * @throws SQLException
     */
    protected String[] lineToArray(ResultSet rs) throws SQLException {
        int colCount = rs.getMetaData().getColumnCount();
        String[] tokens = new String[colCount];
        for (int cc = 0; cc < colCount; cc++) {
            //SQL indexes from 1
            tokens[cc] = rs.getString(cc + 1);
        }
        return tokens;
    }

}
