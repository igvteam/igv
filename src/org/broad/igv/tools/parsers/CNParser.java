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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.tools.parsers;

import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackType;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.readers.AsciiLineReader;

import java.io.IOException;
import java.util.Set;

/**
 * @author jrobinso
 */
public class CNParser extends AbstractParser {

    private static int PROBE_COL = 0;
    int skipColumns;
    // State variables.  This is a serial type parser,  these variables are used to hold temporary
    // state.
    private String chr;
    String lastChr = "";
    int lastPosition = 0;
    private double minValue;
    private double maxValue;
    ResourceLocator resourceLocator;
    Set<String> unsortedChromosomes;
    private int chrColumn;
    private int startColumn;
    private int endColumn;
    private int probeColumn;
    private int firstDataColumn;
    private boolean hasEndLocations;
    private FileType type;
    Genome genome;
    private boolean hasCalls;

    enum FileType {

        IGV, XCN, SNP, CN
    }

    ;

    public CNParser(String file, DataConsumer dataConsumer, Genome genome) {
        this(new ResourceLocator(file), dataConsumer, genome);
    }

    /**
     *
     */
    public CNParser(ResourceLocator locator, DataConsumer dataConsumer, Genome genome) {
        super(dataConsumer);
        this.resourceLocator = locator;
        this.genome = genome;

        String tmp = locator.getPath().toLowerCase();
        tmp = tmp.endsWith(".txt") ? tmp.substring(0, tmp.length() - 4) : tmp;

        hasCalls = tmp.endsWith(".xcn") || tmp.endsWith(".snp");

        if (tmp.endsWith(".igv")) {
            chrColumn = 0;
            startColumn = 1;
            endColumn = 2;
            probeColumn = 3;
            firstDataColumn = 4;
            hasEndLocations = true;
            hasCalls = false;
            type = FileType.IGV;
            setTrackType(TrackType.OTHER);
        } else {
            probeColumn = 0;
            chrColumn = 1;
            startColumn = 2;
            endColumn = -1;
            firstDataColumn = 3;
            hasEndLocations = false;
            hasCalls = tmp.endsWith(".xcn") || tmp.endsWith(".snp");
            if (tmp.endsWith(".cn")) {
                type = FileType.CN;
            } else if (tmp.endsWith(".xcn")) {
                type = FileType.XCN;
            } else {
                type = FileType.SNP;
            }
        }
        skipColumns = hasCalls ? 2 : 1;

    }

    public void parseHeader(String[] tokens) {
        int nDataColumns = (tokens.length - firstDataColumn) / skipColumns;
        String[] headings = new String[nDataColumns];
        for (int i = firstDataColumn; i < tokens.length; i += skipColumns) {
            int idx = (i - firstDataColumn) / skipColumns;
            headings[idx] = tokens[i];
        }
        setHeadings(headings);
    }

    /**
     * @return
     */
    public void parse() throws IOException {

        AsciiLineReader reader = null;
        try {

            lastPosition = 0;

            reader = ParsingUtils.openAsciiReader(resourceLocator);
            String nextLine = null;


            // Infer datatype from extension.  This can be overriden in the
            // comment section
            if (isCopyNumberFileExt(resourceLocator.getPath())) {
                setTrackType(TrackType.COPY_NUMBER);
            } else if (isLOHFileExt(resourceLocator.getPath())) {
                setTrackType(TrackType.LOH);
            }

            // Parse comments, if any
            nextLine = reader.readLine();
            while (nextLine.startsWith("#") || (nextLine.trim().length() == 0)) {
                if (nextLine.length() > 0) {
                    parseComment(nextLine);
                }
                nextLine = reader.readLine();
            }
            parseHeader(nextLine.trim().split("\t"));

            setTrackParameters();

            float[] dataArray = new float[getHeadings().length];


            while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
                String[] tokens = Globals.singleTabMultiSpacePattern.split(nextLine);
                int nTokens = tokens.length;
                if (nTokens == 0) {
                    continue;
                }

                try {

                    chr = (genome == null ? tokens[chrColumn] : genome.getChromosomeAlias(tokens[chrColumn]));
                    if (!chr.equals(lastChr)) {
                        newChromosome();
                    }
                    lastChr = chr;

                    int startPosition = ParsingUtils.parseInt(tokens[startColumn].trim());
                    if (startPosition < lastPosition) {
                        throw new RuntimeException("Error: unsorted file.  .cn files must be sorted by genomic position.");
                    }
                    lastPosition = startPosition;

                    int endPosition = hasEndLocations ? ParsingUtils.parseInt(tokens[endColumn].trim()) : startPosition + 1;

                    // TODO -- compare nTokens with expected number
                    for (int i = firstDataColumn; i < nTokens; i += skipColumns) {
                        int idx = (i - firstDataColumn) / skipColumns;
                        try {
                            dataArray[idx] = Float.parseFloat(tokens[i].trim());
                        } catch (NumberFormatException numberFormatException) {
                            dataArray[idx] = Float.NaN;
                        }
                    }

                    String probe = tokens[probeColumn];

                    getDataConsumer().addData(chr, startPosition, endPosition, dataArray, probe);

                } catch (NumberFormatException e) {
                    System.out.println("Cannot parse: " + nextLine);
                }

            }

            parsingComplete();

        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            if (reader != null) {
                reader.close();
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
        lastPosition = -1;
    }

    private void parsingComplete() {
        getDataConsumer().parsingComplete();
    }

    private boolean isCopyNumberFileExt(String filename) {
        String tmp = (filename.endsWith(".txt") || filename.endsWith(".tab") || filename.endsWith(".xls")
                ? filename.substring(0, filename.length() - 4) : filename);
        return tmp.endsWith(".cn") || tmp.endsWith(".xcn") || tmp.endsWith(".snp");
    }

    private boolean isLOHFileExt(String filename) {
        String tmp = (filename.endsWith(".txt") || filename.endsWith(".tab") || filename.endsWith(".xls")
                ? filename.substring(0, filename.length() - 4) : filename);
        return tmp.endsWith(".loh");
    }
}
