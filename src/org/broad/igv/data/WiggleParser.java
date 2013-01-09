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
package org.broad.igv.data;

import org.apache.commons.math.stat.StatUtils;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.exceptions.ParserException;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.TrackType;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.collections.DownsampledDoubleArrayList;
import org.broad.igv.util.collections.FloatArrayList;
import org.broad.igv.util.collections.IntArrayList;
import org.broad.tribble.readers.AsciiLineReader;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Parser for wiggle and "wiggle-like" formats.
 *
 * @author Enter your name here...
 * @version Enter version here..., 08/10/31
 */
public class WiggleParser {

    private static Logger log = Logger.getLogger(WiggleParser.class);
    private int chrColumn = 0;
    private int startColumn = 1;
    private int endColumn = 2;
    private int dataColumn = 3;

    private enum Type {

        FIXED, VARIABLE, BED_GRAPH, CPG, EXPR
    }

    Genome genome;

    WiggleDataset dataset;
    /**
     * The type of wiggle locator (see UCSC documentation).
     */
    private Type type = Type.BED_GRAPH;

    // State variables.  This is a serial type parser,  these variables are used to hold temporary
    // state.
    private String chr;
    String lastChr = "";
    int lastPosition = 0;
    private int start;
    private int step = 1;
    private int windowSpan = 1;
    private int startBase = 1;   // <- set to zero for zero based coordinates
    IntArrayList startLocations = null;
    IntArrayList endLocations = null;
    FloatArrayList data = null;
    ResourceLocator resourceLocator;
    Set<String> unsortedChromosomes;
    int estArraySize;
    Map<String, Integer> longestFeatureMap = new HashMap();
    // Used to estimate percentiles
    private static final int maxSamples = 1000;
    DownsampledDoubleArrayList sampledData = new DownsampledDoubleArrayList(maxSamples, maxSamples);


    public WiggleParser(ResourceLocator locator) {
        this(locator, null);
    }

    public WiggleParser(ResourceLocator locator, Genome genome) {
        this.genome = genome;
        this.resourceLocator = locator;
        this.estArraySize = estArraySize(locator, genome);
        dataset = new WiggleDataset(genome, locator.getTrackName());

        if (locator.getPath().endsWith("CpG.txt")) {
            type = Type.CPG;
        } else if (locator.getPath().toLowerCase().endsWith(".expr")) {
            //gene_id	bundle_id	chr	left	right	FPKM	FPKM_conf_lo	FPKM_conf_hi
            type = Type.EXPR;
            chrColumn = 2;
            startColumn = 3;
            endColumn = 4;
            dataColumn = 5;
            startBase = 1;
            dataset.setType(TrackType.EXPR);
        }
    }

    private int estArraySize(ResourceLocator locator, Genome genome) {

        int estLines = 100000;
        if (locator.getServerURL() == null) {
            estLines = ParsingUtils.estimateLineCount(locator.getPath());
        }
        int nChromosomes = genome == null ? 24 : genome.getAllChromosomeNames().size();
        return Math.max(1000, (int) (estLines / nChromosomes));

    }

    /**
     * Utility method.  Returns true if this looks like a wiggle locator.  The criteria is to scan
     * the first 100 lines looking for a valid "track" line.  According to UCSC documentation
     * track lines must contain a type attribute,  which must be equal to "wiggle_0".
     *
     * @param file
     * @return
     */
    public static boolean isWiggle(ResourceLocator file) {

        if (file.getPath().endsWith("CpG.txt") || file.getPath().endsWith(".expr") || file.getPath().endsWith(".wig")) {
            return true;
        } else {
            return false;

        }
    }

    public WiggleDataset parse() {

        lastPosition = -1;
        unsortedChromosomes = new HashSet();

        AsciiLineReader reader = null;
        String nextLine = null;
        int lineNumber = 0;

        try {
            reader = ParsingUtils.openAsciiReader(resourceLocator);

            if (type == Type.EXPR) {
                reader.readLine(); // Skip header line
            }

            int position = -1;

            while ((nextLine = reader.readLine()) != null) {
                lineNumber++;

                if (nextLine.startsWith("#") || nextLine.startsWith("data") || nextLine.startsWith("browser") || nextLine.trim().length() == 0) {
                    continue;
                    // Skip
                }


                if (nextLine.startsWith("track") && type != Type.CPG) {
                    type = Type.BED_GRAPH;
                    ParsingUtils.parseTrackLine(nextLine, dataset.getTrackProperties());
                    if (dataset.getTrackProperties().getBaseCoord() == TrackProperties.BaseCoord.ZERO) {
                        this.startBase = 0;
                    }

                } else if (nextLine.startsWith("fixedStep")) {
                    type = Type.FIXED;
                    parseStepLine(nextLine);
                    position = start;
                    if (start < lastPosition) {
                        unsortedChromosomes.add(chr);
                    }

                } else if (nextLine.startsWith("variableStep")) {
                    type = Type.VARIABLE;
                    parseStepLine(nextLine);
                    if (start < lastPosition) {
                        unsortedChromosomes.add(chr);
                    }

                } else {
                    // Must be data
                    String[] tokens = Globals.singleTabMultiSpacePattern.split(nextLine);
                    int nTokens = tokens.length;
                    if (nTokens == 0) {
                        continue;
                    }
                    try {
                        if (type.equals(Type.CPG)) {

                            if (nTokens > 3) {
                                chr = tokens[1].trim();
                                if (!chr.equals(lastChr)) {
                                    changedChromosome(dataset, lastChr);

                                }
                                lastChr = chr;

                                int endPosition = -1;
                                try {
                                    endPosition = Integer.parseInt(tokens[2].trim());
                                } catch (NumberFormatException numberFormatException) {
                                    log.error("Column 2  is not a number");

                                    throw new ParserException("Column 2 must be numeric." + " Found: " + tokens[1],
                                            lineNumber, nextLine);
                                }
                                int startPosition = endPosition - 1;

                                if (startPosition < lastPosition) {
                                    unsortedChromosomes.add(chr);
                                }
                                lastPosition = startPosition;

                                startLocations.add(startPosition);
                                endLocations.add(endPosition);

                                float value = Float.parseFloat(tokens[4].trim());
                                if (tokens[3].trim().equals("R")) {
                                    value = -value;
                                }

                                data.add(value);
                            }
                        } else if (type.equals(Type.BED_GRAPH) || type.equals(Type.EXPR)) {

                            if (nTokens > 3) {
                                chr = tokens[chrColumn].trim();
                                if (!chr.equals(lastChr)) {
                                    changedChromosome(dataset, lastChr);

                                }
                                lastChr = chr;

                                int startPosition = -1;
                                try {
                                    startPosition = Integer.parseInt(tokens[startColumn].trim());
                                } catch (NumberFormatException numberFormatException) {
                                    log.error("Column " + (startColumn + 1) + "  is not a number");

                                    throw new ParserException("Column (startColumn + 1) must be numeric." + " Found: " +
                                            tokens[startColumn],
                                            lineNumber, nextLine);
                                }

                                if (startPosition < lastPosition) {
                                    unsortedChromosomes.add(chr);
                                }
                                lastPosition = startPosition;

                                startLocations.add(startPosition);


                                try {
                                    int endPosition = Integer.parseInt(tokens[endColumn].trim());
                                    endLocations.add(endPosition);
                                    int length = endPosition - startPosition;
                                    updateLongestFeature(length);
                                } catch (NumberFormatException numberFormatException) {
                                    log.error("Column " + (endColumn + 1) + " is not a number");

                                    throw new ParserException("Column " + (endColumn + 1) +
                                            " must be numeric." + " Found: " + tokens[endColumn],
                                            lineNumber, nextLine);
                                }

                                data.add(Float.parseFloat(tokens[dataColumn].trim()));
                            }
                        } else if (type.equals(Type.VARIABLE)) {
                            if (nTokens > 1) {

                                // Per UCSC specification variable and fixed step coordinates are "1" based.
                                // We need to subtract 1 to convert to the internal "zero" based coordinates.
                                int startPosition = Integer.parseInt(tokens[0]) - 1;
                                if (startPosition < lastPosition) {
                                    unsortedChromosomes.add(chr);
                                }
                                lastPosition = startPosition;

                                int end = startPosition + windowSpan;
                                startLocations.add(startPosition);
                                endLocations.add(end);
                                data.add(Float.parseFloat(tokens[1]));
                            }
                        } else {    // Fixed step -- sorting is checked when step line is parsed
                            if (position >= 0) {
                                startLocations.add(position);
                                endLocations.add(position + windowSpan);
                                data.add(Float.parseFloat(tokens[0]));
                            }
                            position += step;
                            lastPosition = position;
                        }

                    } catch (NumberFormatException e) {
                        log.error(e);
                        throw new ParserException(e.getMessage(), lineNumber, nextLine);
                    }


                }

            }

            // The last chromosome
            changedChromosome(dataset, lastChr);

        } catch (ParserException pe) {
            throw (pe);
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


        dataset.sort(unsortedChromosomes);
        dataset.setLongestFeatureMap(longestFeatureMap);

        double[] sd = sampledData.toArray();
        double percent10 = StatUtils.percentile(sd, 10.0);
        double percent90 = StatUtils.percentile(sd, 90.0);
        dataset.setPercent10((float) percent10);
        dataset.setPercent90((float) percent90);

        return dataset;
    }

    private void updateLongestFeature(int length) {
        if (longestFeatureMap.containsKey(chr)) {
            longestFeatureMap.put(chr, Math.max(longestFeatureMap.get(chr), length));
        } else {
            longestFeatureMap.put(chr, length);
        }
    }

    // fixedStep chrom=chrM strt=1 step=1

    private void parseStepLine(String header) {
        String[] tokens = header.split("\\s+");
        for (String token : tokens) {
            String[] keyValue = token.split("=");
            if (keyValue.length >= 2) {
                if (keyValue[0].equalsIgnoreCase("chrom")) {
                    chr = keyValue[1];
                    if (!chr.equals(lastChr)) {
                        changedChromosome(dataset, lastChr);

                    }
                    lastChr = chr;

                } else if (keyValue[0].equalsIgnoreCase("start")) {
                    // Per UCSC specification variable and fixed step coordinates are "1" based.
                    // We need to subtract 1 to convert to the internal "zero" based coordinates.

                    start = Integer.parseInt(keyValue[1]) - startBase;
                    if (start < lastPosition) {
                        unsortedChromosomes.add(chr);
                    }

                } else if (keyValue[0].equalsIgnoreCase("step")) {
                    step = Integer.parseInt(keyValue[1]);
                } else if (keyValue[0].equalsIgnoreCase("span")) {
                    windowSpan = Integer.parseInt(keyValue[1]);
                    updateLongestFeature(windowSpan);
                }

            }
        }
    }


    private void changedChromosome(WiggleDataset dataset, String lastChr) {

        if (startLocations != null && startLocations.size() > 0) {

            String convertedChr = genome == null ? lastChr : genome.getChromosomeAlias(lastChr);
            dataset.addDataChunk(convertedChr, startLocations, endLocations, data);
            //sz = startLocations.size();

            float[] f = data.toArray();
            for (int i = 0; i < f.length; i++) {
                sampledData.add(f[i]);
            }


        }
        startLocations = new IntArrayList(estArraySize);
        endLocations = new IntArrayList(estArraySize);
        data = new FloatArrayList(estArraySize);
        lastPosition = -1;
    }
}
