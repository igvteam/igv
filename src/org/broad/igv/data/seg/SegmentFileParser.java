/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */
package org.broad.igv.data.seg;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.exceptions.ParserException;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.TrackType;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;

/**
 * @author jrobinso
 */
public class SegmentFileParser implements SegFileParser {

    enum Type {
        SEG, BIRDSUITE, NEXUS
    }

    ;

    private static Logger log = Logger.getLogger(SegmentFileParser.class);

    boolean birdsuite = false;
    int sampleColumn = 0;
    int chrColumn = 1;
    int startColumn = 2;
    int endColumn = 3;
    int snpCountColumn = 4;    // Default value
    int dataColumn = 5;        // Default value
    ResourceLocator locator;

    /**
     * Constructs ...
     *
     * @param locator
     */
    public SegmentFileParser(ResourceLocator locator) {
        this.locator = locator;
        if (locator.getPath().toLowerCase().endsWith("birdseye_canary_calls")) {
            birdsuite = true;
        }
    }


    /**
     * Return a map of trackId -> segment datasource
     *
     * @return
     */
    public void loadSegments(SegmentedAsciiDataSet dataset) {

        if (birdsuite) {
            dataset.setTrackType(TrackType.CNV);
        }

        AsciiLineReader reader = null;
        String nextLine = null;
        try {
            reader = ParsingUtils.openAsciiReader(locator);

            // Parse comments, if any
            nextLine = reader.readLine();
            while (nextLine.startsWith("#") || (nextLine.trim().length() == 0)) {
                if (nextLine.length() > 0) {
                    parseComment(nextLine, dataset);
                }
                nextLine = reader.readLine();
            }

            // Read column headings
            String[] headings = nextLine.split("\t");


            if (birdsuite) {
                //sample	sample_index	copy_number	chr	start	end	confidence
                sampleColumn = 0;
                chrColumn = 3;
                startColumn = 4;
                endColumn = 5;
                snpCountColumn = 6;
                dataColumn = 2;
            } else {
                // The data value is always the last column
                dataColumn = headings.length - 1;
                // We assume the snp count is next to last, but not 0-3
                snpCountColumn = dataColumn - 1;
                if (snpCountColumn < 4) {
                    snpCountColumn = -1;
                }
            }

            String[] tokens = new String[headings.length];

            Genome genome = GenomeManager.getInstance().getCurrentGenome();

            while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {

                int nTokens = ParsingUtils.split(nextLine, tokens, '\t');
                if (nTokens > 4) {

                    int start;
                    int end;
                    try {
                        start = Integer.parseInt(tokens[startColumn].trim());
                    }
                    catch (NumberFormatException numberFormatException) {
                        throw new ParserException("Column " + (startColumn + 1) + " must contain a numeric value.",
                                reader.getCurrentLineNumber(), nextLine);
                    }
                    try {
                        end = Integer.parseInt(tokens[endColumn].trim());
                    }
                    catch (NumberFormatException numberFormatException) {
                        throw new ParserException("Column " + (endColumn + 1) + " must contain a numeric value.",
                                reader.getCurrentLineNumber(), nextLine);
                    }

                    String chr = tokens[chrColumn].trim();
                    if (genome != null) {
                        chr = genome.getChromosomeAlias(chr);
                    }


                    String trackId = new String(tokens[sampleColumn].trim());

                    int snpCount = 0;
                    if (snpCountColumn > 0) {
                        try {
                            snpCount = Integer.parseInt(tokens[snpCountColumn]);
                        } catch (NumberFormatException numberFormatException) {

                            // This is an expected condition, nothing needs done.
                        }
                    }

                    try {
                        float value = Float.parseFloat(tokens[dataColumn]);
                        dataset.addSegment(trackId, chr, start, end, value, snpCount);
                    } catch (NumberFormatException numberFormatException) {
                        log.info("Skipping line: " + nextLine);
                    }
                }
            }

        }
        catch (DataLoadException pe) {
            throw pe;
        }
        catch (ParserException pe) {
            throw pe;
        }
        catch (Exception e) {
            if (nextLine != null && reader.getCurrentLineNumber() != 0) {
                throw new ParserException(e.getMessage(), e, reader.getCurrentLineNumber(), nextLine);
            } else {
                throw new RuntimeException(e);
            }
        } finally {
            if (reader == null) {
                reader.close();
            }
        }
    }


    /**
     * Note:  This is an exact copy of the method in GCTDatasetParser.  Refactor to merge these
     * two parsers, or share a common base class.
     *
     * @param comment
     * @param dataset
     */
    private void parseComment(String comment, SegmentedAsciiDataSet dataset) {

        String tmp = comment.substring(1, comment.length());
        String[] tokens = tmp.split("=");
        if (tokens.length != 2) {
            return;
        }
        String key = tokens[0].trim().toLowerCase();
        if (key.equals("type")) {

            try {
                dataset.setTrackType(TrackType.valueOf(tokens[1].trim().toUpperCase()));
            } catch (Exception exception) {

                // Ignore

            }
        }
    }
}
