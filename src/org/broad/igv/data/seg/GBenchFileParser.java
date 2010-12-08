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
import org.broad.igv.synteny.SyntenyMapping;
import org.broad.igv.synteny.SyntenyUtils;
import org.broad.igv.track.TrackType;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;

import java.util.List;
import java.util.Map;

/**
 * @author Enter your name here...
 * @version Enter version here..., 09/01/09
 */
public class GBenchFileParser implements SegFileParser {

    private static Logger log = Logger.getLogger(GBenchFileParser.class);
    String testMappings = "test/data/hg18_to_mm8.regions";

    int chrColumn = 1;
    int startColumn = 3;
    int endColumn = 4;
    int snpCountColumn = 5;
    int ampColumn = 6;
    int delColumn = 7;
    ResourceLocator locator;
    Map<String, List<SyntenyMapping>> mappings;

    /**
     * Constructs ...
     *
     * @param locator
     */
    public GBenchFileParser(ResourceLocator locator) {
        this.locator = locator;


    }

    /**
     * Method description
     *
     * @param dataset
     */
    public void loadSegmentTracks(SegmentedAsciiDataSet dataset) {
        loadSegments(dataset);

    }

    /**
     * Return a map of trackId -> segment datasource
     *
     * @return
     */
    public void loadSegments(SegmentedAsciiDataSet dataset) {
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
            String[] headings = null;

            while ((nextLine = reader.readLine()) != null) {
                if (nextLine.startsWith("AberrationNo")) {
                    headings = nextLine.split("\t");
                    break;
                }
            }
            if (headings == null) {
                throw new DataLoadException("Error: Cannot find data section.", locator.toString());
            }


            Genome genome = GenomeManager.getInstance().getCurrentGenome();

            //
            //if (genome.getId().equals("hg18")) {
            //    mappings = SyntenyUtils.loadMappings(testMappings, true);
            //}

            String trackId = "";
            while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {

                String[] tokens = nextLine.split("\t");
                int nTokens = tokens.length;
                if (nTokens == 1) {
                    trackId = tokens[0];
                } else if (nTokens > 4) {
                    String chr = (genome == null ? tokens[chrColumn].trim() :
                            genome.getChromosomeAlias(tokens[chrColumn].trim()));
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

                    int snpCount = 0;
                    if (snpCountColumn > 0) {
                        try {
                            snpCount = Integer.parseInt(tokens[snpCountColumn]);
                        } catch (NumberFormatException numberFormatException) {

                            // This is an expected condition, nothing needs done.
                        }
                    }


                    float value = Float.NaN;
                    try {
                        value = Float.parseFloat(tokens[ampColumn]);
                        if (value == 0) {
                            value = Float.parseFloat(tokens[delColumn]);
                        }

                    } catch (NumberFormatException numberFormatException) {


                    }

                    if (mappings != null) {
                        List<SyntenyMapping> tmp = mappings.get(chr);
                        if (tmp == null) {
                            log.info("No mappings for chr: " + chr);
                            continue;
                        }
                        List<SyntenyMapping> overlappingMappings = SyntenyUtils.getMappingsOverlapping(tmp, start, end);
                        if(overlappingMappings == null) {
                            log.info("No mapping for: " + chr + ":" + start + "-" + end);
                            continue;
                        }
                        if (overlappingMappings.size() == 1) {
                            SyntenyMapping m1 = overlappingMappings.get(0);
                            int p1 = (int) m1.mapPosition(start);
                            int p2 = (int) m1.mapPosition(end);
                            start = Math.min(p1, p2);
                            end = Math.max(p1, p2);
                            dataset.addSegment(trackId, m1.getToChr(), start, end, value, snpCount);
                        } else {
                            SyntenyMapping m1 = overlappingMappings.get(0);
                            String chr1 = m1.getToChr();
                            int s1 = (int) m1.mapPosition(start);
                            int e1 = m1.getToEnd();
                            if (!m1.getDirection()) {
                                s1 = m1.getToStart();
                                e1 = (int) m1.mapPosition(start);
                            }
                            dataset.addSegment(trackId, chr1, s1, e1, value, snpCount);

                            SyntenyMapping m2 = overlappingMappings.get(overlappingMappings.size() - 1);
                            String chr2 = m2.getToChr();
                            int s2 = (int) m1.mapPosition(end);
                            int e2 = m1.getToEnd();
                            if (!m1.getDirection()) {
                                s2 = m1.getToStart();
                                e2 = (int) m1.mapPosition(end);
                            }
                            dataset.addSegment(trackId, chr2, s2, e2, value, snpCount);

                            for (int i = 1; i < overlappingMappings.size(); i++) {
                                SyntenyMapping m = overlappingMappings.get(i);
                                dataset.addSegment(trackId, m.getToChr(), m.getToStart(), m.getToEnd(), value, snpCount);
                            }
                        }

                    } else {
                        dataset.addSegment(trackId, chr, start, end, value, snpCount);
                    }
                }
            }


        }
        catch (Exception e) {
            e.printStackTrace();
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