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
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.track.FeatureSource;
import org.broad.igv.util.StringUtils;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Types;
import java.util.*;


/**
 * Class for reading data from SQL database, where the lines
 * are of a format we can read with existing codecs
 *
 * @author Jacob Silterra
 * @date 29 May 2012
 */
public class SQLCodecSource extends DBQueryReader<Feature> implements FeatureSource {

    private static Logger log = Logger.getLogger(SQLCodecSource.class);

    public static String UCSC_CHROMO_COL = "chrom";
    public static String UCSC_START_COL = "txStart";
    public static String UCSC_END_COL = "txEnd";

    protected AsciiFeatureCodec codec;

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

    SQLCodecSource(DBProfile.DBTable table, AsciiFeatureCodec codec) {
        super(table);
        this.codec = codec;
        this.binColName = table.getBinColName();
        this.chromoColName = table.getChromoColName();
        this.posStartColName = table.getPosStartColName();
        this.posEndColName = table.getPosEndColName();
        this.startColIndex = table.getStartColIndex();
        this.endColIndex = table.getEndColIndex();

        readHeader();
    }

    /**
     * Read header information from file or database
     * Motivation is mostly for the VCF codec.
     */
    private void readHeader() {

        List<String> headerLines = table.getHeaderLines();
        //If column labels provided in xml spec, use those
        //Otherwise, check the db
        String[] columnLabels;
        if (table.getColumnLabelMap() != null) {
            columnLabels = DBProfile.DBTable.columnMapToArray(table.getColumnLabelMap());
        } else {

            //Preferred method
            /*
            String queryString = String.format("SELECT COLUMN_NAME, ORDINAL_POSITION FROM information_schema.columns WHERE TABLE_NAME = '%s' ORDER BY ORDINAL_POSITION", table.getTableName());

            ResultSet rs = executeQuery(queryString);
            List<String> colLabs = new ArrayList<String>();
            try{
                while(rs.next()){
                    String lab = rs.getString(0);
                    colLabs.add(lab);
                }
            }catch(SQLException e){
                log.error("Error reading column labels", e);
                columnLabels = null;
            } */

            //SQLite doesn't seem to have the information_schema table,
            //and the relevant pragma doesn't have any backwards compatibility guarantees
            //We just query for nothing and read off the column names
            String queryString = String.format("SELECT * FROM %s WHERE 0 = 1", table.getName());
            ResultSet rs = executeQuery(queryString);
            try {
                columnLabels = DBManager.lineToArray(rs, table.getStartColIndex(), table.getEndColIndex(), true);
            } catch (SQLException e) {
                log.error("Error reading column labels", e);
                columnLabels = null;
            }
            DBManager.closeAll(rs);
        }

        if (columnLabels != null) {
            if (headerLines == null) headerLines = new ArrayList<String>(1);
            String columnLine = StringUtils.join(columnLabels, "\t");
            columnLine = "#" + columnLine;
            headerLines.add(columnLine);
        }

        if (headerLines != null) {
            String lines = StringUtils.join(headerLines, "\n");
            byte[] bytes = lines.getBytes();
            InputStream is = new ByteArrayInputStream(bytes);
            LineIterator reader = new LineIteratorImpl(new AsciiLineReader(is));
            try {
                codec.readHeader(reader);
            } catch (IOException e) {
                log.error(e.getMessage(), e);
            }
        }

    }


    /**
     * @param table
     * @return a SQLCodecSource, or null if no appropriate codec found
     */
    public static SQLCodecSource getFromTable(DBProfile.DBTable table) {
        FeatureCodec codec = CodecFactory.getCodec("." + table.getFormat(), GenomeManager.getInstance().getCurrentGenome());
        if (codec != null && codec instanceof AsciiFeatureCodec) {
            return new SQLCodecSource(table, (AsciiFeatureCodec) codec);
        }
        return null;
    }

    public static SQLCodecSource getFromProfile(String profilePath, String tableName) {
        DBProfile dbProfile = DBProfile.parseProfile(profilePath);

        SQLCodecSource source = null;
        for (DBProfile.DBTable table : dbProfile.getTableList()) {
            if (table.getName().equals(tableName)) {
                source = SQLCodecSource.getFromTable(table);
                break;
            }
        }
        return source;
    }

    private String rowToStringLine(ResultSet rs) throws SQLException {
        String[] tokens = rowToStringArray(rs);
        return StringUtils.join(tokens, "\t");
    }

    //TODO We already know how to parse strings, so just turn everything to strings
    //TODO See IParser for better, type-safe way of handling different data sources
    private String[] rowToStringArray(ResultSet rs) throws SQLException {

        String[] tokens;
        if (table.getColumnLabelMap() != null) {
            tokens = DBManager.lineToArray(rs, table.getColumnLabelMap());
        } else {
            tokens = DBManager.lineToArray(rs, startColIndex, endColIndex, false);
        }
        return tokens;
    }

    @Override
    protected Feature processResult(ResultSet rs) throws SQLException {
        String line = rowToStringLine(rs);
        return codec.decode(line);
    }

    /**
     * @param useBinning Whether to query using bin column, for efficiency
     * @throws IOException
     */
    private PreparedStatement generateQueryStatement(boolean useBinning) throws IOException {
        PreparedStatement queryStatement;

        String prependWord = baseQueryString.contains("WHERE") ? " AND " : " WHERE ";
        String queryString = baseQueryString + prependWord + String.format("%s = ? AND ( (%s >= ? AND %s < ?)",
                chromoColName, posStartColName, posStartColName);
        //Don't always have an end position, just assume locations are single base
        if (posEndColName != null) {
            queryString += String.format(" OR (%s < ? AND %s >= ?)", posStartColName, posEndColName);
        }
        queryString += " )";

        String orderClause = "ORDER BY " + posStartColName;

        try {
            if (useBinning) {
                String[] qs = new String[MAX_BINS];
                Arrays.fill(qs, "?");
                String binnedQueryString = queryString + String.format(" AND %s IN (%s) %s", binColName, StringUtils.join(qs, ","), orderClause);
                queryStatement = DBManager.getConnection(locator).prepareStatement(binnedQueryString);
            } else {
                queryStatement = DBManager.getConnection(locator).prepareStatement(queryString + " " + orderClause);
            }
            return queryStatement;
        } catch (SQLException e) {
            log.error("Error initializing query statement", e);
            throw new IOException(e);
        }

    }

    private Iterator query(String chr, int start, int end) throws IOException {

        Set<Integer> bins = null;
        boolean useBinning = false;
        if (binColName != null) {
            bins = calculateBins(start, end);
            useBinning = bins.size() < MAX_BINS;
        }
        PreparedStatement statement = generateQueryStatement(useBinning);

        try {
            statement.clearParameters();
            statement.setString(1, chr);
            statement.setInt(3, end);

            int[] startCols = new int[]{2};
            if (this.posEndColName != null) {
                startCols = new int[]{2, 4, 5};
            }

            for (Integer cc : startCols) {
                statement.setInt(cc, start);
            }

            if (useBinning) {
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
            log.error(e.getMessage(), e);
            throw new IOException(e);
        }

        return loadIterator(statement);
    }


    private static final int SMALLEST_BIN_SIZE = 128 * 1024;

    private Set<Integer> calculateBins(int start, int end) {
        int length = end - start;
        int sweepLength = SMALLEST_BIN_SIZE;
        Set<Integer> bins = new HashSet<Integer>(2 * (length) / SMALLEST_BIN_SIZE);

        while (sweepLength < BINRANGE_MAXEND_512M) {
            int tstStart = Math.max(start - sweepLength / 2, 0);
            while (tstStart < end) {
                bins.add(binFromRange(tstStart, tstStart += sweepLength));
                if (tstStart < 0) {
                    throw new IllegalArgumentException("Overflow while calculating bins");
                }
            }
            sweepLength *= 2;
        }
        bins.add(binFromRange(start, end));

        return bins;
    }

    private static final int BINRANGE_MAXEND_512M = 512 * 1024 * 1024;
    private static final int _binOffsetOldToExtended = 4681;

    /**
     * From http://genomewiki.ucsc.edu/index.php/Bin_indexing_system
     */
    public static int binFromRange(int start, int end) {

        if (start < 0 || end < 0) {
            throw new IllegalArgumentException("start " + start + ", end " + end + " must be > 0");
        }

        boolean extended = false;
        if (end > BINRANGE_MAXEND_512M) {
            extended = true;
        }

        final int binOffsetsExtended[] = {
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
                bin = binOffset + startBin;
                if (extended) {
                    bin += _binOffsetOldToExtended;
                }
                break;
            }
            startBin >>= _binNextShift;
            endBin >>= _binNextShift;
        }

        return bin;
    }


    Iterator iterator() throws IOException {
        String queryString = String.format("%s ORDER BY %s LIMIT %s", baseQueryString, posStartColName, featureWindowSize);
        return loadIterator(executeQuery(queryString));
    }

    public List<String> getSequenceNames() {
        String queryString = String.format("SELECT DISTINCT %s FROM %s", chromoColName, this.getTableName());

        ResultSet results = executeQuery(queryString);
        List<String> names = new ArrayList<String>();
        try {
            while (results.next()) {
                names.add(results.getString(1));
            }
            return names;
        } catch (SQLException e) {
            log.error(e.getMessage(), e);
            throw new RuntimeException(e);
        } finally {
            DBManager.closeAll(results);
        }

    }

    @Override
    public Iterator getFeatures(String chr, int start, int end) throws IOException {
        if (end - start > featureWindowSize) {
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
