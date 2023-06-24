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

package org.broad.igv.gwas;

import org.broad.igv.Globals;
import org.broad.igv.logging.*;
import org.broad.igv.exceptions.ParserException;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;
import java.util.regex.Pattern;

import static java.lang.Math.log10;

/**
 * Parses GWAS PLINK result files
 *
 * @author paananen
 */
public class GWASParser {

    private static final Logger log = LogManager.getLogger(GWASParser.class);

    private ResourceLocator locator;

    Genome genome;
    private GWASColumns columns;

    public Pattern delimiter;

    private int warningCount;

    private static int MAX_WARNING = 20;

    public static boolean isGWASFile(String typeString) {
        return typeString.endsWith("logistic") || typeString.endsWith("linear") || typeString.endsWith("assoc") ||
                typeString.endsWith("qassoc") || typeString.endsWith("gwas");
    }

    public GWASParser(ResourceLocator locator, Genome genome) {
        this.locator = locator;
        this.genome = genome;
        this.columns = new GWASColumns();
        this.warningCount = 0;
    }


    public String[] getColumnHeaders() {
        return this.columns.columnHeaders;
    }


    public GWASData parse() throws IOException {

        BufferedReader reader = null;
        String nextLine = null;
        int rowCounter = 0;
        GWASData data = new GWASData();

        try {
            reader = ParsingUtils.openBufferedReader(locator);

            String headerLine = reader.readLine();

            // Try to determine delimiter pattern -- default is tab, but some older files might use whitespace
            this.delimiter = headerLine.indexOf('\t') > 0 ? Globals.tabPattern : Globals.whitespacePattern;

            if (!this.columns.parseHeader(headerLine, this.delimiter))
                throw new ParserException("Error while parsing header line.", 0, nextLine);

            while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
                nextLine = nextLine.trim();
                rowCounter++;
                GWASFeature f = parseLine(nextLine, rowCounter);
                if (f != null) {
                    data.addFeature(f);
                }
            }
            data.finish();
            return data;

        } catch (Exception e) {
            if (nextLine != null && rowCounter != 0) {
                throw new ParserException(e.getMessage(), e, rowCounter, nextLine);
            } else {
                throw new RuntimeException(e);
            }
        } finally {
            if (reader != null) reader.close();
        }
    }

    /**
     * Parse data from the given text line to {@code GWASData} instance provided
     *
     * @param nextLine
     * @param lineNumber
     * @return Data container, with relevant info
     * @throws ParserException If there is an error parsing the line
     */
    private GWASFeature parseLine(String nextLine, long lineNumber) {

        String[] tokens = this.delimiter.split(nextLine);

        if (tokens.length >= 3) {


            final String posString = tokens[this.columns.locationCol].trim();
            if (posString.indexOf(";") > 0 || posString.length() == 0 || posString.indexOf('x') > 0) {
                if (warningCount < MAX_WARNING) {
                    log.warn(locator.getFileName() + " line number: " + lineNumber + ".  expected numeric position at column " + this.columns.locationCol + " Found " + tokens[this.columns.locationCol]);
                } else if (warningCount == MAX_WARNING) {
                    log.warn("Max warning count excedeed for " + locator.getPath());
                }
                warningCount++;
                return null;
            }

            String chr = tokens[this.columns.chrCol].trim();
            if (genome != null) {
                chr = genome.getCanonicalChrName(chr);
            }

            int position;
            try {
                position = Integer.parseInt(posString);
            } catch (NumberFormatException e) {
                if (warningCount < MAX_WARNING) {
                    log.warn(locator.getFileName() + " line number: " + lineNumber + ".  expected numeric position at column " + this.columns.locationCol + " Found " + tokens[this.columns.locationCol]);
                } else if (warningCount == MAX_WARNING) {
                    log.warn("Max warning count excedeed for " + locator.getPath());
                }
                warningCount++;
                return null;
            }

            // Check if the p-value is NA
            if (!tokens[this.columns.pCol].trim().equalsIgnoreCase("NA")) {
                double p = 0;
                try {
                    // Check for extreme low values
                    String pvalString = tokens[this.columns.pCol];
                    int idx = pvalString.indexOf("E");
                    if (idx > 0) {
                        int exp = Integer.parseInt(pvalString.substring(idx + 1));
                        if (exp < log10(Double.MIN_VALUE)) {
                            p = -1 * exp;
                        }
                    }

                    if (p == 0) {
                        p = Double.parseDouble(tokens[this.columns.pCol]);
                        if (p <= 0) {
                            throw new NumberFormatException();
                        }
                        // Transform to -log10
                        p = -log10(p);
                    }
                } catch (NumberFormatException e) {
                    log.warn("Error parsing line number " + lineNumber + ". Column " + this.columns.pCol + " must be a positive numeric value. Found " + tokens[this.columns.pCol]);
                    return null;
                }
                return new GWASFeature(chr, position, p, nextLine);
            }
        }
        return null;
    }

    /**
     * Stores numerical indexes of relevant columns
     */
    public static class GWASColumns {
        public int locationCol = 2;
        public int chrCol = 1;
        public int pCol = 3;
        public int SNPCol = 0;
        private String[] columnHeaders;

        public boolean hasAllFields() {
            return (this.locationCol >= 0 || this.chrCol >= 0 || this.pCol >= 0 || this.SNPCol >= 0);
        }

        /**
         * Parse a columns string. Based on tokenized columns, populate column numbers to indicate which columns hold chromosome, location, p-value and SNP identifier information.
         *
         * @param headerString
         * @return
         */
        public boolean parseHeader(String headerString, Pattern delimiter) {


            try {
                headerString = headerString.trim();
                String[] headers = delimiter.split(headerString);
                this.columnHeaders = headers;
                int headersSize = headers.length;

                if (headersSize < 4) return false;

                for (int colCounter = 0; colCounter < headersSize; colCounter++) {
                    String header = headers[colCounter];
                    header = header.toLowerCase();

                    // Chromosome column
                    if (header.equals("chr") || header.equals("chromosome") || header.equals("chr_id"))
                        this.chrCol = colCounter;

                    // Nucleotide position column
                    if (header.equals("bp") || header.equals("pos") || header.equals("position") || header.equals("chr_pos"))
                        this.locationCol = colCounter;

                    // p-value column
                    if (header.equals("p") || header.equals("pval") || header.equals("p-value") || header.equals("pvalue") || header.equals("p.value"))
                        this.pCol = colCounter;

                    // SNP identifier column
                    if (header.equals("snp") || header.equals("snps") || header.equals("rs") || header.equals("rsid") || header.equals("rsnum") || header.equals("id") || header.equals("marker") || header.equals("markername"))
                        this.SNPCol = colCounter;
                }
            } catch (Exception e) {
                e.printStackTrace();
            }

            return this.hasAllFields();
        }
    }

}