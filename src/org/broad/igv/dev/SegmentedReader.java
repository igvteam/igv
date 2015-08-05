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

package org.broad.igv.dev;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.data.seg.SegmentedAsciiDataSet;
import org.broad.igv.dev.db.DBManager;
import org.broad.igv.dev.db.DBProfile;
import org.broad.igv.dev.db.WholeTableDBReader;
import org.broad.igv.exceptions.ParserException;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackType;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import htsjdk.tribble.readers.AsciiLineReader;

import java.io.IOException;
import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * Experimental class to explore using the same class as a reader for
 * both a SQL DB and a file.
 *
 * @author Jacob Silterra
 * @date 2012 Aug 30
 */
public class SegmentedReader {

    private static Logger log = Logger.getLogger(SegmentedReader.class);
    private Genome genome;
    private SegmentedAsciiDataSet dataset;
    private String[] headings;

    private boolean birdsuite = false;
    private int sampleColumn = 0;
    private int chrColumn = 1;
    private int startColumn = 2;
    private int endColumn = 3;
    //int snpCountColumn = 4;    // Default value
    private int dataColumn = 5;
    private ResourceLocator locator;

    enum Type {
        SEG, BIRDSUITE, NEXUS
    }

    private class SegmentedDBReader extends WholeTableDBReader<SegmentedAsciiDataSet> {

        private SegmentedDBReader(DBProfile.DBTable table) {
            super(table);
        }

        @Override
        protected SegmentedAsciiDataSet processResultSet(ResultSet rs) throws SQLException {

            IParser<ResultSet, Integer> parser = IParserFactory.getIndexParser(rs);
            readHeader(rs);

            while (rs.next()) {
                parseLine(parser, rs);
            }
            DBManager.closeAll(rs);

            dataset.sortLists();
            return dataset;
        }

        private int readHeader(ResultSet rs) throws SQLException {
            dataset = new SegmentedAsciiDataSet(genome);
            /*sampleColumn = 0;
            chrColumn = 1;
            startColumn = 2;
            endColumn = 3;*/
            dataColumn = rs.getMetaData().getColumnCount() - 1;

            headings = new String[rs.getMetaData().getColumnCount()];
            for (int cc = 0; cc < rs.getMetaData().getColumnCount(); cc++) {
                headings[cc] = rs.getMetaData().getColumnLabel(cc + 1);
            }

            return 0;
        }

    }

    public SegmentedReader(ResourceLocator locator) {
        this(locator, null);
    }

    public SegmentedReader(ResourceLocator locator, Genome genome) {
        this.locator = locator;
        this.genome = genome;

        if (locator.getPath().toLowerCase().endsWith("birdseye_canary_calls")) {
            birdsuite = true;
        }
    }

    /**
     * Load segmented data from a file could be remote)
     * Return a map of trackId -> segment datasource
     *
     * @return
     */
    public SegmentedAsciiDataSet loadFromFile() {

        dataset = new SegmentedAsciiDataSet(genome);

        if (birdsuite) {
            dataset.setTrackType(TrackType.CNV);
        }

        AsciiLineReader reader = null;
        String nextLine = null;
        int lineNumber = 0;

        IParser<String[], Integer> parser = IParserFactory.getIndexParser(new String[0]);

        try {
            reader = ParsingUtils.openAsciiReader(locator);
            lineNumber = readHeader(reader);

            while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
                lineNumber++;
                String[] tokens = Globals.tabPattern.split(nextLine, -1);
                parseLine(parser, tokens);
            }

        } catch (IOException e) {
            if (nextLine != null && lineNumber != 0) {
                throw new ParserException(e.getMessage(), e, lineNumber, nextLine);
            } else {
                throw new RuntimeException(e);
            }
        } finally {
            if (reader != null) {
                reader.close();
            }
        }

        dataset.sortLists();
        return dataset;
    }


    public SegmentedAsciiDataSet loadFromDB(DBProfile.DBTable table) {
        SegmentedDBReader reader = new SegmentedDBReader(table);
        return reader.load();
    }

    /**
     * Note:  This is an exact copy of the method in ExpressionFileParser.  Refactor to merge these
     * two parsers, or share a common base class.
     *
     * @param comment
     * @param dataset
     */
    private void parseComment(String comment, SegmentedAsciiDataSet dataset) {

        String tmp = comment.substring(1, comment.length());
        if (tmp.startsWith("track")) {
            ParsingUtils.parseTrackLine(tmp, dataset.getTrackProperties());
        } else {
            String[] tokens = tmp.split("=");

            String key = tokens[0].trim().toLowerCase();
            if (key.equals("type")) {
                if (tokens.length != 2) {
                    return;
                }
                try {
                    dataset.setTrackType(TrackType.valueOf(tokens[1].trim().toUpperCase()));
                } catch (Exception exception) {
                    log.error("Unknown track type: " + tokens[1].trim().toUpperCase());
                }
            }

        }
    }

    private int readHeader(AsciiLineReader reader) throws IOException {
        // Parse comments, if any
        String nextLine = reader.readLine();
        int lineNumber = 0;
        while (nextLine.startsWith("#") || (nextLine.trim().length() == 0)) {
            lineNumber++;
            if (nextLine.length() > 0) {
                parseComment(nextLine, dataset);
            }
            nextLine = reader.readLine();
        }

        // Read column headings
        headings = nextLine.split("\t");

        if (birdsuite) {
            //sample	sample_index	copy_number	chr	start	end	confidence
            sampleColumn = 0;
            dataColumn = 2;
            chrColumn = 3;
            startColumn = 4;
            endColumn = 5;

        } else {
            sampleColumn = 0;
            chrColumn = 1;
            startColumn = 2;
            endColumn = 3;
            dataColumn = headings.length - 1;
        }

        return lineNumber;
    }

    private <TContainer> void parseLine(IParser<TContainer, Integer> parser, TContainer rs) throws ParserException {

        int nTokens = parser.size(rs);
        if (nTokens > 4) {

            String sample = parser.getString(rs, sampleColumn);
            String chr = parser.getString(rs, chrColumn);
            int start = parser.getInt(rs, startColumn);
            int end = parser.getInt(rs, endColumn);
            float value = parser.getFloat(rs, dataColumn);


            if (genome != null) {
                chr = genome.getChromosomeAlias(chr);
            }

            String trackId = sample;

            StringBuffer desc = null;
            if (birdsuite) {
                desc = new StringBuffer();
                desc.append("<br>");
                desc.append(headings[6]);
                desc.append("=");
                desc.append(parser.getString(rs, 6));
            } else {
                if (nTokens > 4) {
                    desc = new StringBuffer();
                    for (int i = 4; i < parser.size(rs) - 1; i++) {
                        desc.append("<br>");
                        desc.append(headings[i]);
                        desc.append(": ");
                        desc.append(parser.getString(rs, i));
                    }
                }
            }
            String description = desc == null ? null : desc.toString();
            dataset.addSegment(trackId, chr, start, end, value, description);
        }

    }

}
