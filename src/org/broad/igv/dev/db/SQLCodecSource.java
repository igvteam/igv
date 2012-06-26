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
import org.broad.tribble.AsciiFeatureCodec;
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
import java.sql.*;
import java.util.*;


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
    public static String UCSC_START_COL = "txStart";
    public static String UCSC_END_COL = "txEnd";


    protected AsciiFeatureCodec codec;

    protected PreparedStatement queryStatement;
    protected PreparedStatement binnedQueryStatement;

    /**
     * The name of the column with chromosome names
     * Default is correct value for UCSC genes
     */
    protected String chromoColName = UCSC_CHROMO_COL;
    /**
     * The name of the column of positions that we query over.
     * Default is correct value for UCSC genes
     */
    protected String posStartColName = UCSC_START_COL;


    /**
     * Name of the column which contains
     * the end location
     */
    protected String posEndColName = UCSC_END_COL;

    /**
     * Some databases use a "bin" column to speed up queries.
     * See http://genomewiki.ucsc.edu/index.php/Bin_indexing_system
     * or doi: 10.1101/gr.229102 Genome Res. 2002. 12: 996-1006
     */
    protected String binColName;

    /**
     * We start reading data from this column, by default is 1 (SQL indexes from 1)
     * If there are leading index columns but the table structure is otherwise
     * similar to a file format, can skip by setting startColIndex to
     * the first data column.
     */
    protected int startColIndex = 1;

    /**
     * Maximum column number to read, by default is Integer.MAX_VALUE
     * Columns beyond this are ignored. Inclusive.
     */
    protected int endColIndex = Integer.MAX_VALUE;
    private int featureWindowSize = (int) 1e6;

    private static final int MAX_BINS = 20;

    public SQLCodecSource(ResourceLocator locator, AsciiFeatureCodec codec, String table) {
        super(locator, table);
        this.codec = codec;
    }

    public SQLCodecSource(ResourceLocator locator, AsciiFeatureCodec codec, String table,
                          String chromoColName, String posStartColName, String posEndColName, int startColIndex, int endColIndex) {
        this(locator, codec, table);
        this.chromoColName = chromoColName;
        this.posStartColName = posStartColName;
        this.posEndColName = posEndColName;
        this.startColIndex = startColIndex;
        this.endColIndex = endColIndex;
    }

    /**
     * Retrieve a reader from the XML profile located at {@code profilePath}.
     * TODO If {@code tableName == null}, the user is prompted to choose a table from the list
     *
     * @param profilePath
     * @param tableName
     * @return
     */
    public static List<SQLCodecSource> getFromProfile(String profilePath, String tableName) {
        ResourceLocator dbLocator = DBManager.getStoredConnection(profilePath);
        InputStream profileStream = null;
        try {
            profileStream = new FileInputStream(profilePath);
            Document document = Utilities.createDOMDocumentFromXmlStream(profileStream);
            NodeList nodes = document.getElementsByTagName("table");
            List<SQLCodecSource> sources = new ArrayList<SQLCodecSource>(nodes.getLength());

            for (int tnum = 0; tnum < nodes.getLength(); tnum++) {
                Node n = nodes.item(tnum);
                NamedNodeMap attr = n.getAttributes();
                String tabName = attr.getNamedItem("name").getTextContent();
                if (tableName == null || tableName.equals(tabName)) {

                    String chromoColName = attr.getNamedItem("chromoColName").getTextContent();
                    String posStartColName = attr.getNamedItem("posStartColName").getTextContent();
                    String posEndColName = attr.getNamedItem("posEndColName").getTextContent();
                    String format = attr.getNamedItem("format").getTextContent();
                    String startColString = Utilities.getNullSafe(attr, "startColIndex");
                    String endColString = Utilities.getNullSafe(attr, "endColIndex");
                    String binColName = Utilities.getNullSafe(attr, "binColName");
                    int startColIndex = startColString != null ? Integer.parseInt(startColString) : 1;
                    int endColIndex = endColString != null ? Integer.parseInt(endColString) : Integer.MAX_VALUE;
                    AsciiFeatureCodec codec = CodecFactory.getCodec("." + format, GenomeManager.getInstance().getCurrentGenome());
                    SQLCodecSource source = new SQLCodecSource(dbLocator, codec, tabName, chromoColName, posStartColName, posEndColName, startColIndex, endColIndex);
                    source.binColName = binColName;
                    sources.add(source);
                }
            }

            return sources;

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
        int colCount = Math.min(rs.getMetaData().getColumnCount(), endColIndex) - startColIndex + 1;
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
        //Include feature iff = (feature.start >= start AND feature.start < end)
        //OR (feature.start < start AND feature.end >= start);
        String queryString = String.format("%s WHERE %s = ? AND ( (%s >= ? AND %s < ?) OR (%s < ? AND %s >= ?))",
                baseQueryString, chromoColName, posStartColName, posStartColName, posStartColName, posEndColName);

        try {
            queryStatement = DBManager.getConnection(locator).prepareStatement(queryString);

            if(binColName != null){
                String[] qs = new String[MAX_BINS];
                Arrays.fill(qs, "?");
                String binnedQueryString = queryString + String.format(" AND %s IN (%s)", binColName, StringUtils.join(qs, ','));
                binnedQueryStatement = DBManager.getConnection(locator).prepareStatement(binnedQueryString);
            }
        } catch (SQLException e) {
            log.error(e);
            throw new IOException(e);
        }

    }

    private CloseableTribbleIterator query(String chr, int start, int end) throws IOException {
        initQueryStatement();
        PreparedStatement statement = queryStatement;
        Set<Integer> bins = calculateBins(start, end);
        //System.out.println("number of bins: " + bins.size());
        if(bins.size() < MAX_BINS && binnedQueryStatement != null){
            statement = binnedQueryStatement;
        }

        try {
            statement.clearParameters();
            statement.setString(1, chr);
            statement.setInt(3, end);
            int[] cols = new int[]{2, 4, 5};
            for(Integer cc: cols){
                statement.setInt(cc, start);
            }

            if (statement == binnedQueryStatement) {
                int qnum = 6;
                for (Integer bin : bins) {
                    statement.setInt(qnum, bin);
                    qnum++;
                }

                for (; qnum <= statement.getParameterMetaData().getParameterCount(); qnum++) {
                    statement.setNull(qnum, Types.INTEGER);
                }
            }

        } catch (SQLException e) {
            log.error(e);
            throw new IOException(e);
        }

        return loadIterator(statement);
    }


    private static final int SMALLEST_BIN_SIZE = 128*1024;

    private Set<Integer> calculateBins(int start, int end) {
        int length = end - start;
        int sweepLength = SMALLEST_BIN_SIZE;
        Set<Integer> bins = new HashSet<Integer>(2*(length) / SMALLEST_BIN_SIZE);

        while(sweepLength < BINRANGE_MAXEND_512M){
            int tstStart = Math.max(start - sweepLength/2, 0);
            while(tstStart < end){
                bins.add(binFromRange(tstStart, tstStart += sweepLength));
            }
            sweepLength *= 2;
        }
        bins.add(binFromRange(start, end));

        return bins;
    }

    private static final int BINRANGE_MAXEND_512M = 512*1024*1024;
    private static final int _binOffsetOldToExtended = 4681;

    /**
     * From http://genomewiki.ucsc.edu/index.php/Bin_indexing_system
     */
    public static int binFromRange(int start, int end) {

        if(start < 0 || end < 0){
            throw new IllegalArgumentException("start "+start+", end "+ end + " must be > 0");
        }

        boolean extended = false;
        if(end > BINRANGE_MAXEND_512M){
            extended = true;
        }

        final int binOffsetsExtended[] ={
               4096 + 512 + 64 + 8 + 1,
                      512 + 64 + 8 + 1,
                            64 + 8 + 1,
                                 8 + 1,
                                     1,
                                     0};

        int[] binOffsets = Arrays.copyOfRange(binOffsetsExtended,
                extended ? 0 : 1, binOffsetsExtended.length);


        /** How much to shift to get to first bin. */
        final int _binFirstShift = 17;
        /** How much to shift to get to next larger bin.*/
        final int _binNextShift = 3;

        int startBin = start;
        int endBin = end - 1;

        startBin >>= _binFirstShift;
        endBin >>= _binFirstShift;
        int bin = -1;

        for (int binOffset : binOffsets) {

            if (startBin == endBin) {
                bin= binOffset + startBin;
                if(extended){
                    bin += _binOffsetOldToExtended;
                }
                break;
            }
            startBin >>= _binNextShift;
            endBin >>= _binNextShift;
        }

        return bin;
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
