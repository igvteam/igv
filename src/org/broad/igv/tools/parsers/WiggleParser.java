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
package org.broad.igv.tools.parsers;

//~--- non-JDK imports --------------------------------------------------------

import com.iontorrent.utils.StringTools;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.TrackType;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.readers.AsciiLineReader;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

/**
 * Class description
 *
 * @author Enter your name here...
 * @version Enter version here..., 08/10/31
 */
public class WiggleParser {

    private enum Type {

        fixedStep, variableStep, bed, cpg
    }

    private DataConsumer dataConsumer;
    /**
     * The type of wiggle locator (see UCSC documentation).
     */
    private Type type = Type.bed;
    // State variables.  This is a serial type parser,  these variables are used to hold temporary
    // state.
    String trackLine = null;
    String nextLine = null;
    private String chr;
    String lastChr = "";
    int lastPosition = 0;
    private int start = 0;
    private int step = 1;
    private int span = 1;    // <- ignored for now
    private int startBase = 1;   // <- set to zero for zero based coordinates
    private double minValue;
    private double maxValue;
    ResourceLocator resourceLocator;
    Set<String> unsortedChromosomes;
    Genome genome;

    /**
     *
     */
    public WiggleParser(String file, DataConsumer dataConsumer, Genome genome) {
        this.resourceLocator = new ResourceLocator(file);
        this.dataConsumer = dataConsumer;
        this.genome = genome;

        parseHeader();

        String[] trackNames = {resourceLocator.getTrackName()};

        // TODO -- total hack to get Manuel's file parsed quickly.  Revisit (obviously);
        if (resourceLocator.getPath().endsWith(".ewig") || resourceLocator.getPath().endsWith(".ewig.gz")
                || resourceLocator.getPath().endsWith("ewig.map")) {
            trackNames = new String[5];
            trackNames[4] = resourceLocator.getTrackName();
            trackNames[0] = "A";
            trackNames[1] = "C";
            trackNames[2] = "G";
            trackNames[3] = "T";
        }

        dataConsumer.setTrackParameters(TrackType.OTHER, trackLine, trackNames);


        // Parse track line, if any, to get the coordinate convention
        if (trackLine != null) {
            TrackProperties props = new TrackProperties();
            ParsingUtils.parseTrackLine(trackLine, props);
            TrackProperties.BaseCoord convention = props.getBaseCoord();
            if (convention == TrackProperties.BaseCoord.ZERO) {
                startBase = 0;
            }
        }


        if (resourceLocator.getPath().endsWith("CpG.txt")) {
            type = Type.cpg;
        }

    }

    /**
     * @return the dataConsumer
     */
    public DataConsumer getDataConsumer() {
        return dataConsumer;
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
        AsciiLineReader reader = null;
        try {
            reader = ParsingUtils.openAsciiReader(file);
            String nextLine = null;
            int lineNo = 0;
            while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
                if (nextLine.startsWith("track") && nextLine.contains("wiggle_0")) {
                    return true;
                }
                if (lineNo++ > 100) {
                    break;
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
            return false;
        } finally {
            if (reader != null) {
                reader.close();
            }
        }
        return false;
    }

    private void parseHeader() {
        AsciiLineReader reader = null;

        // The DataConsumer interface takes an array of data per position, however wig
        // files contain a single data point.  Create an "array" once that can
        // be resused

        try {

            reader = ParsingUtils.openAsciiReader(resourceLocator);

            while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {

                // Skip comment lines
                if (nextLine.startsWith("#") || nextLine.startsWith("data") || nextLine.startsWith(
                        "browser") || nextLine.trim().length() == 0) {
                    continue;
                }

                if (nextLine.startsWith("track")) {
                    trackLine = nextLine;

                } else {
                    return;
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (reader != null) {
                reader.close();
            }
        }
    }

    /**
     * @return
     */
    public void parse() throws IOException {

        lastPosition = 0;
        unsortedChromosomes = new HashSet();

        if (resourceLocator.getPath().endsWith("ewig.map")) {
            startBase = 0;
            span = 1;
            type = Type.variableStep;
            String parent = new File(resourceLocator.getPath()).getParent();
            Map<String, String> fileMap = parseEwigList(resourceLocator.getPath());
            for (Map.Entry<String, String> entry : fileMap.entrySet()) {
                chr = entry.getKey();
                File f = new File(parent, entry.getValue());
                parseFile(f.getAbsolutePath());
            }


        } else {
            parseFile(resourceLocator.getPath());
        }
        parsingComplete();

    }

    private Map<String, String> parseEwigList(String path) throws IOException {
        BufferedReader reader = null;
        try {
            reader = ParsingUtils.openBufferedReader(path);
            String nextLine;
            LinkedHashMap<String, String> fileMap = new LinkedHashMap();
            while ((nextLine = reader.readLine()) != null) {
                String[] tokens = Globals.whitespacePattern.split(nextLine);
                fileMap.put(tokens[0], tokens[1]);
            }
            return fileMap;

        } finally {
            if (reader != null) reader.close();
        }
    }

    private void parseFile(String path) throws IOException {
        AsciiLineReader reader = null;

        // The DataConsumer interface takes an array of data per position, however wig
        // files contain a single data point.  Create an "array" once that can
        // be resused
        float[] dataArray = null;

        try {

            reader = ParsingUtils.openAsciiReader(new ResourceLocator(path));

            int position = -1;
            while ((nextLine = reader.readLine()) != null) {

                // Skip comment lines
                if (nextLine.startsWith("#") || nextLine.startsWith("data") || nextLine.startsWith("browser") || nextLine.trim().length() == 0) {
                    continue;
                }


                if (nextLine.startsWith("track")) {
                    // BED by default
                    type = Type.bed;
                    dataConsumer.setType("bed");
                    //DatasetParserUtils.parseTrackLine(nextLine, dataset.getTrackProperties());

                } else if (nextLine.startsWith("fixedStep")) {
                    type = Type.fixedStep;
                    dataConsumer.setType("fixedStep");

                    parseStepLine(nextLine);
                    position = start;
                    if (start < lastPosition) {
                        unsortedChromosomes.add(chr);
                    }


                } else if (nextLine.startsWith("variableStep")) {
                    type = Type.variableStep;
                    dataConsumer.setType("variableStep");

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
                        if (type.equals(Type.cpg)) {
                            if (nTokens > 3) {

                                // TODO -- validation on the data array length (# of data columns).
                                if (dataArray == null) {
                                    dataArray = new float[nTokens - 3];
                                }

                                chr = tokens[1].trim();
                                if (!chr.equals(lastChr)) {
                                    newChromosome();
                                }
                                lastChr = chr;

                                int endPosition = -1;
                                try {
                                    endPosition = Integer.parseInt(tokens[2].trim());
                                } catch (NumberFormatException numberFormatException) {
                                    System.err.println("Column 2  is not a number");

                                    throw new RuntimeException("Column 2 must be numeric." + " Found: " + tokens[1]);
                                }
                                int startPosition = endPosition - 1;

                                if (startPosition < lastPosition) {
                                    unsortedChromosomes.add(chr);
                                }
                                lastPosition = startPosition;


                                float value = Float.parseFloat(tokens[4].trim());
                                if (tokens[3].trim().equals("R")) {
                                    value = -value;
                                }
                                dataArray[0] = value;
                                getDataConsumer().addData(chr, startPosition, endPosition, dataArray, null);
                            }
                        } else if (type.equals(Type.bed)) {
                            if (nTokens > 3) {

                                // TODO -- validation on the data array length (# of data columns).
                                if (dataArray == null) {
                                    dataArray = new float[nTokens - 3];
                                }

                                chr = (genome == null ? tokens[0].trim() : genome.getChromosomeAlias(tokens[0].trim()));
                                if (!chr.equals(lastChr)) {
                                    newChromosome();
                                }
                                lastChr = chr;

                                int startPosition = Integer.parseInt(tokens[1].trim());
                                if (startPosition < lastPosition) {
                                    unsortedChromosomes.add(chr);
                                }
                                lastPosition = startPosition;

                                int endPosition = Integer.parseInt(tokens[2].trim());

                                for (int i = 0; i < dataArray.length; i++) {
                                    dataArray[i] = Float.parseFloat(tokens[3 + i].trim());
                                }

                                getDataConsumer().addData(chr, startPosition, endPosition, dataArray, null);
                            }
                        } else if (type.equals(Type.variableStep)) {
                            if (nTokens > 1) {

                                if (dataArray == null) {
                                    dataArray = new float[nTokens - 1];
                                }


                                // Per UCSC specification variable and fixed step coordinates are "1" based.
                                // We need to subtract 1 to convert to the internal "zero" based coordinates.
                                int startPosition = Integer.parseInt(tokens[0]) - startBase;
                                if (startPosition < lastPosition) {
                                    unsortedChromosomes.add(chr);
                                }
                                lastPosition = startPosition;

                                int endPosition = startPosition + span;

                                for (int i = 0; i < dataArray.length; i++) {
                                    dataArray[i] = Float.parseFloat(tokens[1 + i].trim());
                                }

                                getDataConsumer().addData(chr, startPosition, endPosition, dataArray, null);
                            }
                        } else {    // Fixed step -- sorting is checked when step line is parsed
                            if (position >= 0) {

                                if (dataArray == null) {
                                    dataArray = new float[nTokens];
                                }

                                int endPosition = position + span;

                                for (int i = 0; i < dataArray.length; i++) {
                                    dataArray[i] = Float.parseFloat(tokens[i].trim());
                                }

                                getDataConsumer().addData(chr, position, endPosition, dataArray, null);

                            }
                            position += step;
                            lastPosition = position;

                        }

                    } catch (NumberFormatException e) {
                        System.out.println("Cannot parse: " + nextLine);
                    }
                }
            }

        } finally {
            if (reader != null) {
                reader.close();
            }
        }
    }

    // fixedStep chrom=chrM strt=1 step=1

    private void parseStepLine(String header) {
        String[] tokens = header.split("\\s+");
        for (String token : tokens) {
            String[] keyValue = token.split("=");
            if (keyValue.length >= 2) {
                if (keyValue[0].equalsIgnoreCase("chrom")) {
                    chr = (genome == null ? keyValue[1] : genome.getChromosomeAlias(keyValue[1]));
                    if (!chr.equals(lastChr)) {
                        newChromosome();
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
                    span = Integer.parseInt(keyValue[1]);
                }

            }
        }
    }

    /**
     * Method description
     *
     * @return
     */
    public double getMinValue() {
        return minValue;
    }

    /**
     * Method description
     *
     * @return
     */
    public double getMaxValue() {
        return maxValue;
    }

    private void newChromosome() {
        //getDataConsumer().newChromosome(chr);
        //lastPosition = 0;
    }

    private void parsingComplete() {
        getDataConsumer().parsingComplete();
    }
}
