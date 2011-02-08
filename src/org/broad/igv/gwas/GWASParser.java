/*
 * Copyright (c) 2007-2009 by The Broad Institute, Inc. and the Massachusetts Institute of Technology.
 * All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL), Version 2.1 which
 * is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR WARRANTIES OF
 * ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT
 * OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR
 * RESPECTIVE TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES OF
 * ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES, ECONOMIC
 * DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER THE BROAD OR MIT SHALL
 * BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE
 * FOREGOING.
 */
package org.broad.igv.gwas;

import org.apache.log4j.Logger;
import org.broad.igv.exceptions.ParserException;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.readers.AsciiLineReader;

import java.io.FileInputStream;
import java.io.IOException;

import static java.lang.Math.log10;

/**
 * Parses GWAS PLINK result files
 *
 * @author paananen
 */
public class GWASParser {

    private static final Logger log = Logger.getLogger(GWASParser.class);
    private ResourceLocator locator;

    private int locationCol = -1;
    private int chrCol = -1;
    private int pCol = -1;
    private int SNPCol = -1;


    public ResourceLocator getLocator() {
        return locator;
    }

    public void setLocator(ResourceLocator locator) {
        this.locator = locator;
    }


    public GWASParser(ResourceLocator locator) {
        this.locator = locator;

    }

    /**
     * Parse a header string. Based on tokenized header, populate column numbers to indicate which columns hold chromosome, location, p-value and SNP identifier information.
     *
     * @param headerString
     * @return
     */
    public boolean parseHeader(String headerString) {

        boolean parsingSuccessful = true;

        headerString = headerString.trim();


        String[] headers = new String[1000];

        int headersSize = ParsingUtils.splitSpaces(headerString, headers);

        if (headersSize < 4)
            parsingSuccessful = false;

        int colCounter = 0;
        for (int i = 0; i < headersSize; i++) {
            String header = headers[i];
            header = header.toLowerCase();

            // Chromosome column
            if (header.equals("chr") || header.equals("chromosome"))
                this.chrCol = colCounter;

            // Nucleotide position column
            if (header.equals("bp") || header.equals("pos") || header.equals("position"))
                this.locationCol = colCounter;

            // p-value column
            if (header.equals("p") || header.equals("pval") || header.equals("p-value") || header.equals("pvalue") || header.equals("p.value"))
                this.pCol = colCounter;

            // SNP identifier column
            if (header.equals("snp") || header.equals("rs") || header.equals("rsid") || header.equals("rsnum") || header.equals("id") || header.equals("marker") || header.equals("markername"))
                this.SNPCol = colCounter;

            colCounter++;

        }

        if (this.locationCol < 0 || this.chrCol < 0 || this.pCol < 0 || this.SNPCol < 0)
            parsingSuccessful = false;

        return parsingSuccessful;
    }


    /**
     * Parses and populates description cache from a GWAS result file. Cache will be filled with data points surrounding the given query data point.
     *
     * @param gData          GWASData object
     * @param hitChr         Chromosome of the query data point
     * @param hitLocation    Nucleotide location of the query data point
     * @param searchStartRow Result file row where populating the cache will start
     * @return
     * @throws IOException
     */
    public GWASData parseDescriptions(GWASData gData, String hitChr, long hitLocation, int searchStartRow) throws IOException {


        FileInputStream fs = null;

        AsciiLineReader reader = null;
        String nextLine = null;
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        boolean hitFound = false;
        int cacheSize = gData.getDescriptionCache().getMaxSize();

        try {
            fs = new FileInputStream(locator.getPath());
            fs.getChannel().position(0);
            reader = new AsciiLineReader(fs);

            // Parse header line
            String headerLine = reader.readLine();
            if (!parseHeader(headerLine))
                throw new ParserException("Error while parsing header line.", reader.getCurrentLineNumber(), nextLine);

            gData.getDescriptionCache().setHeaderTokens(headerLine);

            int rowCounter = 0;
            // Number of lines to cache after the hit (tries to estimate that half of the cache will be filled with data from data before the query data point, and half with data after the query data point
            int cacheAfter = cacheSize / 2;
            int cacheCounter = 0;

            while (cacheCounter < cacheAfter && (nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {

                nextLine = nextLine.trim();
                rowCounter++;

                if (rowCounter >= searchStartRow) {
                    String[] tokens = new String[1000];
                    ParsingUtils.splitSpaces(nextLine, tokens);

                    if (tokens.length > 1) {

                        // Get chromosome
                        String chr = genome.getChromosomeAlias(tokens[chrCol].trim());

                        // Get nucleotide position
                        int start;
                        try {
                            start = Integer.parseInt(tokens[locationCol].trim());
                        } catch (NumberFormatException e) {
                            throw new ParserException("Column " + locationCol + " must be a numeric value.", reader.getCurrentLineNumber(), nextLine);
                        }

                        // See if chr and nucleotide position match to our query data point
                        if (chr.equals(hitChr) && start == hitLocation) {
                            hitFound = true;

                        }

                        // Add descriptions to cache
                        if (hitFound) {
                            cacheCounter++;
                        }
                        gData.getDescriptionCache().add(chr, start, nextLine);
                    }
                }
            }


        } catch (
                ParserException e
                )

        {
            throw e;
        } catch (
                Exception e
                )

        {
            if (nextLine != null && reader.getCurrentLineNumber() != 0) {
                throw new ParserException(e.getMessage(), e, reader.getCurrentLineNumber(), nextLine);
            } else {
                throw new RuntimeException(e);
            }
        } finally

        {
            reader.close();
            fs.close();
        }

        return gData;

    }

    /**
     * Find a description line from the result file for a single data point. This is more efficiently achieved using the DescriptionCache(), but could be useful for certain types of data sets.
     *
     * @param hitChr         Chromosome of the query data point
     * @param hitLocation    Nucleotide location of the query data point
     * @param searchStartRow Result file row where populating the cache will start
     * @return
     * @throws IOException
     */
    public String parseDescriptionLine(String hitChr, double hitLocation, int searchStartRow) throws IOException {


        FileInputStream fs = null;

        AsciiLineReader reader = null;
        String nextLine = null;
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        boolean hitFound = false;
        String resultLine = "";


        try {
            fs = new FileInputStream(locator.getPath());
            fs.getChannel().position(0);
            reader = new AsciiLineReader(fs);

            // Parse header line
            String headerLine = reader.readLine();
            if (!parseHeader(headerLine))
                throw new ParserException("Error while parsing header line.", reader.getCurrentLineNumber(), nextLine);

            // Tokenize header for creation of tokenized descriptions
            headerLine = headerLine.trim();
            String[] headers = new String[1000];
            int headersSize = ParsingUtils.splitSpaces(headerLine, headers);


            int rowCounter = 0;

            while (!hitFound && (nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {


                nextLine = nextLine.trim();
                rowCounter++;
                if (rowCounter >= searchStartRow) {
                    String[] tokens = new String[100];
                    ParsingUtils.splitSpaces(nextLine, tokens);

                    if (tokens.length > 1) {

                        String chr = genome.getChromosomeAlias(tokens[chrCol].trim());

                        int start;

                        try {
                            start = Integer.parseInt(tokens[locationCol].trim());
                        } catch (NumberFormatException e) {
                            throw new ParserException("Column " + locationCol + " must be a numeric value.", reader.getCurrentLineNumber(), nextLine);
                        }

                        if (chr.equals(hitChr) && start == hitLocation) {
                            hitFound = true;
                            for (int i = 0; i < headersSize; i++) {
                                resultLine += headers[i] + ": " + tokens[i] + "<br>";

                            }
                        }
                    }
                }
            }


            return resultLine;
        } catch (
                ParserException e
                )

        {
            throw e;
        } catch (
                Exception e
                )

        {
            if (nextLine != null && reader.getCurrentLineNumber() != 0) {
                throw new ParserException(e.getMessage(), e, reader.getCurrentLineNumber(), nextLine);
            } else {
                throw new RuntimeException(e);
            }
        } finally

        {
            reader.close();
            fs.close();
        }


    }

    public GWASData parse() throws IOException {


        FileInputStream fs = null;

        AsciiLineReader reader = null;
        String nextLine = null;
        Genome genome = GenomeManager.getInstance().getCurrentGenome();


        try {
            fs = new FileInputStream(locator.getPath());
            fs.getChannel().position(0);
            reader = new AsciiLineReader(fs);


            // Parse header line
            String headerLine = reader.readLine();
            if (!parseHeader(headerLine))
                throw new ParserException("Error while parsing header line.", reader.getCurrentLineNumber(), nextLine);

            GWASData gData = new GWASData();

            int rowCounter = 0;
            int indexCounter = 0;
            int addedValuesCounter = 0;

            while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {


                nextLine = nextLine.trim();
                rowCounter++;

                String[] tokens = new String[100];
                ParsingUtils.splitSpaces(nextLine, tokens);

                if (tokens.length > 1) {

                    //String chr = ParsingUtils.convertChrString(tokens[chrCol].trim());
                    String chr = genome.getChromosomeAlias(tokens[chrCol].trim());

                    int start;

                    try {
                        start = Integer.parseInt(tokens[locationCol].trim());
                    } catch (NumberFormatException e) {
                        throw new ParserException("Column " + locationCol + " must be a numeric value.", reader.getCurrentLineNumber(), nextLine);
                    }

                    // Check if the p-value is NA
                    if (!tokens[pCol].trim().equals("NA")) {
                        float p;

                        try {
                            p = Float.parseFloat(tokens[pCol].trim());
                            // Transform to -log10
                            p = (float) -log10((double) p);


                        } catch (NumberFormatException e) {
                            throw new ParserException("Column " + pCol + " must be a numeric value.", reader.getCurrentLineNumber(), nextLine);
                        }


                        gData.addLocation(chr, start);
                        gData.addValue(chr, p);

                        indexCounter++;
                        addedValuesCounter++;

                        int indexSize = 10000;
                        if (indexCounter == indexSize) {
                            gData.getFileIndex().add((int) reader.getPosition());

                            indexCounter = 0;

                        }


                    }

                }


            }


            return gData;
        } catch (
                ParserException e
                )

        {
            throw e;
        } catch (
                Exception e
                )

        {
            if (nextLine != null && reader.getCurrentLineNumber() != 0) {
                throw new ParserException(e.getMessage(), e, reader.getCurrentLineNumber(), nextLine);
            } else {
                throw new RuntimeException(e);
            }
        } finally

        {
            reader.close();
            fs.close();
        }


    }


}