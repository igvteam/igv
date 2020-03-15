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

import org.apache.log4j.Logger;
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

    private static final Logger log = Logger.getLogger(GWASParser.class);

    private ResourceLocator locator;

    Genome genome;
    private GWASColumns columns;

    public static boolean isGWASFile(String typeString) {
        return typeString.endsWith(".logistic") || typeString.endsWith(".linear") || typeString.endsWith(".assoc") ||
                typeString.endsWith(".qassoc") || typeString.endsWith(".gwas");
    }

    public GWASParser(ResourceLocator locator, Genome genome) {
        this.locator = locator;
        this.genome = genome;
        this.columns = new GWASColumns();
    }


    public String[] getColumnHeaders() {
        return this.columns.columnHeaders;
    }


    public Map<String, List<GWASFeature>> parse() throws IOException {

        BufferedReader reader = null;
        String nextLine = null;
        int rowCounter = 0;
        Map<String, List<GWASFeature>> features = new HashMap<>();

        try {
            reader = ParsingUtils.openBufferedReader(locator);

            String headerLine = reader.readLine();

            if (!this.columns.parseHeader(headerLine))
                throw new ParserException("Error while parsing header line.", 0, nextLine);

            while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
                nextLine = nextLine.trim();
                rowCounter++;
                GWASFeature f = null;
                f = parseLine(nextLine, rowCounter);
                List<GWASFeature> featureList = features.get(f.chr);
                if (featureList == null) {
                    featureList = new ArrayList<>();
                    features.put(f.chr, featureList);
                }
                featureList.add(f);
            }

            // Sort features by position
            for (List<GWASFeature> featureList : features.values()) {
                featureList.sort(Comparator.comparingInt(o -> o.position));
            }

            return features;

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

        String[] tokens = nextLine.split("\\t|( +)");

        if (tokens.length >= 3) {
            String chr = tokens[this.columns.chrCol].trim();
            if (genome != null) {
                chr = genome.getCanonicalChrName(chr);
            }

            int position;
            try {
                position = Integer.parseInt(tokens[this.columns.locationCol].trim());
            } catch (NumberFormatException e) {
                throw new ParserException("Column " + this.columns.locationCol + " must be a numeric value.", lineNumber, nextLine);
            }

            // Check if the p-value is NA
            if (!tokens[this.columns.pCol].trim().equalsIgnoreCase("NA")) {
                double p;
                try {
                    p = Double.parseDouble(tokens[this.columns.pCol]);
                    if (p <= 0) {
                        throw new NumberFormatException();
                    }
                    // Transform to -log10
                    p = -log10(p);
                } catch (NumberFormatException e) {
                    throw new ParserException("Column " + this.columns.pCol + " must be a positive numeric value. Found " + tokens[this.columns.pCol], lineNumber, nextLine);
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
        public boolean parseHeader(String headerString) {

            Pattern whitespacePattern = Pattern.compile("\\s+");
            try {
                headerString = headerString.trim();
                String[] headers = headerString.split("\\t|( +)");
                this.columnHeaders = headers;
                int headersSize = headers.length;

                if (headersSize < 4) return false;

                for (int colCounter = 0; colCounter < headersSize; colCounter++) {
                    String header = headers[colCounter];
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
                }
            } catch (Exception e) {
                e.printStackTrace();
            }

            return this.hasAllFields();
        }
    }

    public static void main(String[] args) throws IOException {
        GWASParser parser = new GWASParser(new ResourceLocator("test/data/gwas/smallp.gwas"), null);
        Map<String, List<GWASFeature>> data = parser.parse();
        for (List<GWASFeature> features : data.values()) {
            for (GWASFeature f : features) {
                double val = f.value;
                System.out.println(val);
            }
        }
    }
}