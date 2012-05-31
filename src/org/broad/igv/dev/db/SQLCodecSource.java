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
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.track.FeatureSource;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.Utilities;
import org.broad.tribble.CloseableTribbleIterator;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.sql.Blob;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Class for reading data from SQL database, where the lines
 * are of a format we can read with existing codecs
 *
 * @author Jacob Silterra
 * @date 29 May 2012
 */
public class SQLCodecSource extends DBReader<Feature> implements FeatureSource {

    private static Logger log = Logger.getLogger(SQLCodecSource.class);

    public static String UCSC_CHROMO_COL = "chrom";
    public static String UCSC_POS_COL = "txStart";

    protected FeatureCodec codec;

    protected PreparedStatement queryStatement;

    /**
     * The name of the column with chromosome names
     * Default is correct value for UCSC genes
     */
    protected String chromoColName = UCSC_CHROMO_COL;
    /**
     * The name of the column of positions that we query over.
     * Default is correct value for UCSC genes
     */
    protected String posColName = UCSC_POS_COL;

    /**
     * We start reading data from this column, by default is 1 (SQL indexes from 1)
     * If there are leading index columns but the table structure is otherwise
     * similar to a file format, can skip by setting startColIndex to
     * the first data column.
     */
    protected int startColIndex = 1;
    private int featureWindowSize = (int) 1e6;

    public SQLCodecSource(ResourceLocator locator, FeatureCodec codec, String table) {
        super(locator, table);
        this.codec = codec;
    }

    public SQLCodecSource(ResourceLocator locator, FeatureCodec codec, String table, String chromoColName, String posColName, int startColIndex) {
        this(locator, codec, table);
        this.chromoColName = chromoColName;
        this.posColName = posColName;
        this.startColIndex = startColIndex;
    }

    /**
     * Retrieve a reader from the XML profile located at {@code profilePath}.
     * TODO If {@code tableName == null}, the user is prompted to choose a table from the list
     *
     * @param profilePath
     * @param tableName
     * @return
     */
    public static SQLCodecSource getFromProfile(String profilePath, String tableName) {
        ResourceLocator dbLocator = DBManager.getStoredConnection(profilePath);
        InputStream profileStream = null;
        try {
            profileStream = new FileInputStream(profilePath);
            Document document = Utilities.createDOMDocumentFromXmlStream(profileStream);
            NodeList nodes = document.getElementsByTagName("table");
            if (nodes.getLength() != 1 && tableName == null) {
                throw new IllegalArgumentException("Found " + nodes.getLength() + " tables but didn't specify a name");
            }

            for (int tnum = 0; tnum < nodes.getLength(); tnum++) {
                Node n = nodes.item(tnum);
                NamedNodeMap attr = n.getAttributes();
                String tabName = attr.getNamedItem("name").getTextContent();
                if (tableName == null || tableName.equals(tabName)) {

                    String chromoColName = attr.getNamedItem("chromoColName").getTextContent();
                    String posColName = attr.getNamedItem("posColName").getTextContent();
                    String format = attr.getNamedItem("format").getTextContent();
                    String startColString = Utilities.getNullSafe(attr, "startColIndex");
                    int startColIndex = startColString != null ? Integer.parseInt(startColString) : 1;
                    FeatureCodec codec = CodecFactory.getCodec("." + format, GenomeManager.getInstance().getCurrentGenome());
                    return new SQLCodecSource(dbLocator, codec, tableName, chromoColName, posColName, startColIndex);
                }
            }

        } catch (Exception e) {
            e.printStackTrace();
            return null;
        } finally {
            try {
                if (profileStream != null) profileStream.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        return null;

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
        int colCount = rs.getMetaData().getColumnCount() - startColIndex + 1;
        String[] tokens = new String[colCount];
        String s;
        int sqlCol;
        for (int cc = 0; cc < colCount; cc++) {

            //SQL indexes from 1
            //Have to parse blobs specially, otherwise we get the pointer as a string
            sqlCol = cc + startColIndex;
            String clazz = rs.getMetaData().getColumnClassName(sqlCol).toLowerCase();
            if (clazz.contains("blob") || clazz.equalsIgnoreCase("[b")) {
                Blob b = rs.getBlob(sqlCol);
                s = new String(b.getBytes(1l, (int) b.length()));
            } else {
                s = rs.getString(sqlCol);
            }
            tokens[cc] = s;
        }
        return tokens;
    }

    /**
     * Create the prepared statement. Idempotent.
     *
     * @throws IOException
     */
    private void initQueryStatement() throws IOException {
        if (queryStatement != null) {
            return;
        }
        String queryString = String.format("%s WHERE %s = ? AND %s >= ? AND %s < ?",
                baseQueryString, chromoColName, posColName, posColName);
        try {
            queryStatement = DBManager.getConnection(locator).prepareStatement(queryString);
        } catch (SQLException e) {
            log.error(e);
            throw new IOException(e);
        }

    }

    private CloseableTribbleIterator query(String chr, int start, int end) throws IOException {
        initQueryStatement();
        try {
            queryStatement.clearParameters();
            queryStatement.setString(1, chr);
            queryStatement.setInt(2, start);
            queryStatement.setInt(3, end);
        } catch (SQLException e) {
            log.error(e);
            throw new IOException(e);
        }
        return loadIterator(queryStatement);
    }

    CloseableTribbleIterator iterator() throws IOException {
        String queryString = String.format("%s LIMIT %s", baseQueryString, featureWindowSize);
        return loadIterator(queryString);
    }

    private void close() throws IOException {
        try {
            queryStatement.close();
            queryStatement = null;
            DBManager.closeConnection(locator);
        } catch (SQLException e) {
            log.error(e);
            throw new IOException(e);
        }
    }

    public List<String> getSequenceNames() {
        String queryString = String.format("SELECT DISTINCT %s FROM %s", chromoColName, table);

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
    public Iterator getFeatures(String chr, int start, int end) throws IOException {
        if (start - end > featureWindowSize) {
            return null;
        }
        return query(chr, start, end);
    }

    @Override
    public List<LocusScore> getCoverageScores(String chr, int start, int end, int zoom) {
        return null; //TODO
    }

    @Override
    public int getFeatureWindowSize() {
        return featureWindowSize;
    }

    @Override
    public void setFeatureWindowSize(int size) {
        this.featureWindowSize = size;
    }
}
