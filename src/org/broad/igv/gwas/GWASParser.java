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
package org.broad.igv.gwas;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.exceptions.ParserException;
import org.broad.igv.feature.genome.Genome;
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
    Genome genome;

    public static boolean isGWASFile(String typeString) {
        return typeString.endsWith(".logistic") || typeString.endsWith(".linear") || typeString.endsWith(".assoc") ||
                typeString.endsWith(".qassoc") || typeString.endsWith(".gwas");
    }


    public GWASParser(ResourceLocator locator, Genome genome) {
        this.locator = locator;
        this.genome = genome;
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


        String[] headers = Globals.singleTabMultiSpacePattern.split(headerString);
        int headersSize = headers.length;

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
     * TODO -- This method is nearly identical to "parse()", and is apparently only used to support the popup text.
     * <p/>
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

        boolean hitFound = false;
        int cacheSize = gData.getDescriptionCache().getMaxSize();
        int rowCounter = 0;

        try {
            fs = new FileInputStream(locator.getPath());
            fs.getChannel().position(0);
            reader = new AsciiLineReader(fs);

            // Parse header line
            String headerLine = reader.readLine();
            if (!parseHeader(headerLine))
                throw new ParserException("Error while parsing header line.", 0, nextLine);

            gData.getDescriptionCache().setHeaderTokens(headerLine);


            // Number of lines to cache after the hit (tries to estimate that half of the cache will be filled with data from data before the query data point, and half with data after the query data point
            int cacheAfter = cacheSize / 2;
            int cacheCounter = 0;

            while (cacheCounter < cacheAfter && (nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {

                nextLine = nextLine.trim();
                rowCounter++;

                if (rowCounter >= searchStartRow) {
                    String[] tokens = Globals.singleTabMultiSpacePattern.split(nextLine);

                    if (tokens.length > 1) {


                        // Check if the p-value is NA
                        if (!tokens[pCol].trim().equals("NA")) {
                            double p;

                            try {
                                p = Double.parseDouble(tokens[pCol].trim());
                                // Transform to -log10
                                p = -log10(p);


                            } catch (NumberFormatException e) {
                                throw new ParserException("Column " + pCol + " must be a numeric value.", rowCounter, nextLine);
                            }


                            // Get chromosome
                            String chr = genome.getChromosomeAlias(tokens[chrCol].trim());

                            // Get nucleotide position
                            int start;
                            try {
                                start = Integer.parseInt(tokens[locationCol].trim());
                            } catch (NumberFormatException e) {
                                throw new ParserException("Column " + locationCol + " must be a numeric value.", rowCounter, nextLine);
                            }

                            // See if chr and nucleotide position match to our query data point
                            if (chr.equals(hitChr) && start == hitLocation) {
                                hitFound = true;

                            }

                            // Add descriptions to cache
                            if (hitFound) {
                                cacheCounter++;
                            }
                            gData.getDescriptionCache().add(chr, start, p, nextLine);
                        }
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
            if (nextLine != null && rowCounter != 0) {
                throw new ParserException(e.getMessage(), e, rowCounter, nextLine);
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

    public GWASData parse() throws IOException {


        FileInputStream fs = null;

        AsciiLineReader reader = null;
        String nextLine = null;
        int rowCounter = 0;

        try {
            fs = new FileInputStream(locator.getPath());
            fs.getChannel().position(0);
            reader = new AsciiLineReader(fs);


            // Parse header line
            String headerLine = reader.readLine();
            if (!parseHeader(headerLine))
                throw new ParserException("Error while parsing header line.", 0, nextLine);

            GWASData gData = new GWASData();

            int indexCounter = 0;
            int addedValuesCounter = 0;

            while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {


                nextLine = nextLine.trim();
                rowCounter++;

                String[] tokens = Globals.singleTabMultiSpacePattern.split(nextLine);

                if (tokens.length > 1) {

                    //String chr = ParsingUtils.convertChrString(tokens[chrCol].trim());
                    String chr = genome.getChromosomeAlias(tokens[chrCol].trim());

                    int start;

                    try {
                        start = Integer.parseInt(tokens[locationCol].trim());
                    } catch (NumberFormatException e) {
                        throw new ParserException("Column " + locationCol + " must be a numeric value.", rowCounter, nextLine);
                    }

                    // Check if the p-value is NA
                    if (!tokens[pCol].trim().equalsIgnoreCase("NA")) {
                        double p;

                        try {
                            p = Double.parseDouble(tokens[pCol]);
                            if (p <= 0) {
                                throw new NumberFormatException();
                            }
                            // Transform to -log10
                            p = -log10(p);

                        } catch (NumberFormatException e) {
                            throw new ParserException("Column " + pCol + " must be a positive numeric value. Found " + tokens[pCol], rowCounter, nextLine);
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

        } catch (Exception e) {
            if (nextLine != null && rowCounter != 0) {
                throw new ParserException(e.getMessage(), e, rowCounter, nextLine);
            } else {
                throw new RuntimeException(e);
            }
        } finally {
            reader.close();
            fs.close();
        }


    }


}