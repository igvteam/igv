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
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.CloseableTribbleIterator;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.FeatureReader;

import java.io.IOException;
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
public class SQLCodecReader extends DBReader<Feature> implements FeatureReader {

    private static Logger log = Logger.getLogger(SQLCodecReader.class);

    public static String UCSC_CHROMO_COL = "chrom";
    public static String UCSC_POS_COL = "txStart";

    protected FeatureCodec codec;


    /**
     * The name of the column with chromosome names
     * Default is correct value for UCSC genes
     */
    protected String chromoCol = UCSC_CHROMO_COL;
    /**
     * The name of the column of positions that we query over.
     * Default is correct value for UCSC genes
     */
    protected String posCol = UCSC_POS_COL;

    public SQLCodecReader(ResourceLocator locator, FeatureCodec codec, String table) {
        super(locator, table);
        this.codec = codec;
    }

    public SQLCodecReader(ResourceLocator locator, FeatureCodec codec, String table, String chromoCol, String posCol) {
        this(locator, codec, table);
        this.chromoCol = chromoCol;
        this.posCol = posCol;
    }

    @Override
    protected Feature processResult(ResultSet rs) throws SQLException {
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

    @Override
    public CloseableTribbleIterator query(String chr, int start, int end) throws IOException {
        String queryString = String.format("%s WHERE %s = '%s' AND %s >= %d AND %s < %d",
                baseQueryString, chromoCol, chr, posCol, start, posCol, end);
        return loadIterator(queryString);
    }

    @Override
    public CloseableTribbleIterator iterator() throws IOException {
        //TODO Should we allow this at all?
        String queryString = String.format("%s LIMIT 100000", baseQueryString);
        return loadIterator(queryString);
    }

    @Override
    public void close() throws IOException {
        DBManager.closeConnection(this.locator);
    }

    @Override
    public List<String> getSequenceNames() {
        String queryString = String.format("SELECT DISTINCT %s FROM %s", chromoCol, table);

        ResultSet results = loadResultSet(queryString);
        List<String> names = new ArrayList<String>();
        try {
            while (results.next()) {
                names.add(results.getString(1));
            }
            return names;
        } catch (SQLException e) {
            log.error(e);
            throw new RuntimeException("Error getting sequence names: " + e);
        }

    }

    @Override
    public Object getHeader() {
        return null;
    }
}
