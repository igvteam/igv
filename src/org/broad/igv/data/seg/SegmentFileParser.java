/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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
import org.broad.igv.track.TrackType;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.readers.AsciiLineReader;

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

            String[] tokens = new String[headings.length];

            while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {

                int nTokens = ParsingUtils.split(nextLine, tokens, '\t');
                if (nTokens > 4) {


                    int start;
                    int end;
                    try {
                        start = (int) Double.parseDouble(tokens[startColumn].trim());
                    } catch (NumberFormatException numberFormatException) {
                        throw new ParserException("Column " + (startColumn + 1) + " must contain a numeric value.",
                                reader.getCurrentLineNumber(), nextLine);
                    }
                    try {
                        end = (int) Double.parseDouble(tokens[endColumn].trim());
                    } catch (NumberFormatException numberFormatException) {
                        throw new ParserException("Column " + (endColumn + 1) + " must contain a numeric value.",
                                reader.getCurrentLineNumber(), nextLine);
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
