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

/*
 * GisticFilesParser.java
 *
 * Created on June 22, 2007, 8:32 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package org.broad.igv.feature;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.exceptions.ParserException;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.GisticTrack;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import htsjdk.tribble.readers.AsciiLineReader;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;

/**
 * @author jrobinso
 */
public class GisticFileParser {

    private static Logger log = Logger.getLogger(GisticFileParser.class);

    /**
     * Method description
     *
     * @param locator
     * @return
     */
    public static GisticTrack loadData(ResourceLocator locator) {

        // TODO -- handle remote resource
        AsciiLineReader reader = null;

        String nextLine = null;
        int rowCounter = 0;
        try {
            Genome genome = GenomeManager.getInstance().getCurrentGenome();

            reader = ParsingUtils.openAsciiReader(locator);

            // Parse header
            reader.readLine();
            rowCounter++;

            // Parse data
            // parameters to track chromosome breaks. Used in wholeGenome case.
            GisticTrack track = new GisticTrack(locator);
            track.setName(locator.getTrackName());
            List<GisticScore> scores = new ArrayList();

            nextLine = reader.readLine();
            rowCounter++;

            while ((nextLine != null) && (nextLine.length() > 0)) {
                String[] tokens = nextLine.split("\t");

                GisticScore.Type type = getType(tokens[0].trim());
                String chr = genome.getChromosomeAlias(tokens[1].trim());

                int start = -1;
                try {
                    start = Integer.parseInt(tokens[2]);
                }
                catch (NumberFormatException ne) {
                    throw new ParserException("Column 3 must be a numeric value.", rowCounter, nextLine);
                }

                int end = -1;
                try {
                    end = Integer.parseInt(tokens[3]);
                }
                catch (NumberFormatException ne) {
                    throw new ParserException("Column 4 must be a numeric value.", rowCounter, nextLine);
                }

                float qValue = 0;
                try {
                    qValue = Float.parseFloat(tokens[4]);
                } catch (NumberFormatException numberFormatException) {
                    if (tokens[4].toLowerCase().equals("inf")) {
                        qValue = Float.POSITIVE_INFINITY;
                    } else {
                        qValue = Float.NaN;
                    }
                }

                try {
                    float gScore = Float.parseFloat(tokens[5]);
                    scores.add(new GisticScore(chr, start, end, qValue, gScore, type));

                }
                catch (Exception ex) {
                    log.error("Skipping line: " + nextLine, ex);
                }

                nextLine = reader.readLine();
                rowCounter++;
            }

            if (scores.isEmpty()) {
                throw new RuntimeException("No gistic scores were found.");
            } else {
                track.setScores(scores);
                computeWholeGenome(track, genome);
                return track;
            }

        } catch (FileNotFoundException e) {
            log.error("File not found: " + locator.getPath());
            throw new RuntimeException(e);
        }
        catch (ParserException e) {
            throw new RuntimeException(e);
        }
        catch (Exception e) {
            log.error("Exception when loading: " + locator.getPath(), e);
            if (rowCounter != 0) {
                throw new ParserException(e.getMessage(), rowCounter, nextLine);
            } else {
                throw new RuntimeException(e);
            }
        } finally {
            if (reader != null) {
                reader.close();
            }
        }

    }

    /**
     * Compute features for "chr All".  This is probably not the best way to do this,
     * definitely not good for large tracks.  Need it now for mutation data.
     *
     * @param track
     */
    public static void computeWholeGenome(GisticTrack track, Genome genome) {
        int unit = 1000;    // KB

        String chrAll = "All";
        List<GisticScore> allFeatures = new ArrayList(1000);
        if (genome == null) {
            throw new RuntimeException("Unknown genome: " + genome.getId());
        }
        genome.getLongChromosomeNames();
        long offset = 0;
        for (String chr : genome.getLongChromosomeNames()) {
            int chrLength = genome.getChromosome(chr).getLength();
            int chrStart = (int) (offset / unit);
            allFeatures.add(new GisticScore(chrAll, chrStart, chrStart, 0, 0, GisticScore.Type.AMP));
            allFeatures.add(new GisticScore(chrAll, chrStart, chrStart, 0, 0, GisticScore.Type.DEL));
            List<GisticScore> ampScores = track.getAmpScores(chr);
            if (ampScores != null) {
                for (GisticScore m : ampScores) {

                    int start = genome.getGenomeCoordinate(chr, m.getStart());
                    int end = genome.getGenomeCoordinate(chr, m.getEnd());

                    allFeatures.add(new GisticScore(chrAll, start, end,
                            (float) m.getQValue(),
                            (float) m.getGScore(),
                            m.getType()));
                }
            }
            List<GisticScore> delScores = track.getDelScores(chr);
            if (delScores != null) {
                for (GisticScore m : delScores) {
                    int start = (int) ((offset + m.getStart()) / unit);
                    int end = (int) ((offset + m.getEnd()) / unit);

                    allFeatures.add(new GisticScore(chrAll, start, end,
                            (float) m.getQValue(),
                            (float) m.getGScore(),
                            m.getType()));
                }
            }

            offset += chrLength;
        }
        track.setScores(allFeatures);

    }

    static GisticScore.Type getType(String typeString) {
        if (typeString.toUpperCase().equals("AMP")) {
            return GisticScore.Type.AMP;
        } else if (typeString.toUpperCase().equals("DEL")) {
            return GisticScore.Type.DEL;
        } else {
            throw new IllegalArgumentException("Unkown score type: " + typeString);
        }
    }
}
