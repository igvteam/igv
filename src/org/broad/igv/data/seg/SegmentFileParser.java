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

package org.broad.igv.data.seg;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.exceptions.ParserException;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackType;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import htsjdk.tribble.readers.AsciiLineReader;

/**
 * Example
 * CCLE_name	chrom	loc.start	loc.end	num.mark	seg.mean
 * A2780_OVARY	1	51598	88465	17	1.0491
 * A2780_OVARY	1	218569	55606710	30512	0.0248
 * A2780_OVARY	1	55606801	61331401	4399	-0.0431
 *
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
    //int snpCountColumn = 4;    // Default value
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
    public SegmentedAsciiDataSet loadSegments(ResourceLocator locator, Genome genome) {

        SegmentedAsciiDataSet dataset = new SegmentedAsciiDataSet(genome);

        if (birdsuite) {
            dataset.setTrackType(TrackType.CNV);
        }

        AsciiLineReader reader = null;
        String nextLine = null;
        int lineNumber = 0;
        try {
            reader = ParsingUtils.openAsciiReader(locator);

            // Parse comments, if any
            nextLine = reader.readLine();
            while (nextLine.startsWith("#") || (nextLine.trim().length() == 0)) {
                lineNumber++;
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

            while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
                lineNumber++;

                String[] tokens = Globals.tabPattern.split(nextLine, -1);
                int nTokens = tokens.length;
                if (nTokens > 4) {
                    int start;
                    int end;
                    try {
                        start = ParsingUtils.parseInt(tokens[startColumn].trim());
                    } catch (NumberFormatException numberFormatException) {
                        throw new ParserException("Column " + (startColumn + 1) + " must contain a numeric value.",
                                lineNumber, nextLine);
                    }
                    try {
                        end = ParsingUtils.parseInt(tokens[endColumn].trim());
                    } catch (NumberFormatException numberFormatException) {
                        throw new ParserException("Column " + (endColumn + 1) + " must contain a numeric value.",
                                lineNumber, nextLine);
                    }

                    String chr = tokens[chrColumn].trim();
                    if (genome != null) {
                        chr = genome.getChromosomeAlias(chr);
                    }


                    String trackId = new String(tokens[sampleColumn].trim());

                    StringBuffer desc = null;
                    if (birdsuite) {
                        desc = new StringBuffer();
                        desc.append("<br>");
                        desc.append(headings[6]);
                        desc.append("=");
                        desc.append(tokens[6]);
                    } else {
                        if (tokens.length > 4) {
                            desc = new StringBuffer();
                            for (int i = 4; i < headings.length - 1; i++) {
                                desc.append("<br>");
                                desc.append(headings[i]);
                                desc.append(": ");
                                desc.append(tokens[i]);
                            }
                        }
                    }


                    try {
                        float value = Float.parseFloat(tokens[dataColumn]);
                        String description = desc == null ? null : desc.toString();
                        dataset.addSegment(trackId, chr, start, end, value, description);
                    } catch (NumberFormatException numberFormatException) {
                        // log.info("Skipping line: " + nextLine);
                    }
                }
            }

        } catch (DataLoadException pe) {
            throw pe;
        } catch (ParserException pe) {
            throw pe;
        } catch (Exception e) {
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
}
