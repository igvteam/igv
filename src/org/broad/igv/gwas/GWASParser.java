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
import org.broad.igv.Globals;
import org.broad.igv.exceptions.ParserException;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.StringUtils;
import htsjdk.tribble.readers.AsciiLineReader;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

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

        AsciiLineReader reader = null;
        String nextLine = null;

        boolean hitFound = false;
        int cacheSize = gData.getDescriptionCache().getMaxSize();
        int rowCounter = 0;

        try {
            reader = ParsingUtils.openAsciiReader(locator);

            // Parse columns line
            String headerLine = reader.readLine();
            if (!this.columns.parseHeader(headerLine))
                throw new ParserException("Error while parsing columns line.", 0, nextLine);

            gData.getDescriptionCache().setHeaderTokens(headerLine);


            // Number of lines to cache after the hit (tries to estimate that half of the cache will be filled with data from data before the query data point, and half with data after the query data point
            int cacheAfter = cacheSize / 2;
            int cacheCounter = 0;

            while (cacheCounter < cacheAfter && (nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {

                nextLine = nextLine.trim();
                rowCounter++;

                if (rowCounter >= searchStartRow) {

                    GWASEntry entry = parseLine(nextLine, rowCounter);
                    if(entry == null) continue;

                    // See if chr and nucleotide position match to our query data point
                    if (entry.chr.equals(hitChr) && entry.start == hitLocation) {
                        hitFound = true;
                    }

                    // Add descriptions to cache
                    if (hitFound) {
                        cacheCounter++;
                    }
                    gData.getDescriptionCache().add(entry.chr, entry.start, entry.p, nextLine);

                }
            }

        } catch (ParserException e) {
            throw e;
        } catch (Exception e){
            if (nextLine != null && rowCounter != 0) {
                throw new ParserException(e.getMessage(), e, rowCounter, nextLine);
            } else {
                throw new RuntimeException(e);
            }
        } finally{
            if(reader != null) reader.close();
        }

        return gData;

    }

    public GWASData parse() throws IOException {

        AsciiLineReader reader = null;
        String nextLine = null;
        int rowCounter = 0;

        Set<String> chromos = new HashSet<String>();
        GWASEntry lastEntry = null;

        try {
            reader = ParsingUtils.openAsciiReader(locator);

            String headerLine = reader.readLine();
            if (!this.columns.parseHeader(headerLine))
                throw new ParserException("Error while parsing columns line.", 0, nextLine);

            GWASData gData = new GWASData();

            int indexCounter = 0;

            while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {

                nextLine = nextLine.trim();
                rowCounter++;

                GWASEntry entry = parseLine(nextLine, rowCounter);
                if (entry == null) continue;

                gData.addLocation(entry.chr, entry.start);
                gData.addValue(entry.chr, entry.p);

                //Check that file is sorted
                if(lastEntry != null){
                    if(entry.chr.equals(lastEntry.chr)){
                        if(entry.start < lastEntry.start){
                            throw new ParserException("File is not sorted, found start position lower than previous", rowCounter);
                        }
                    }else{
                        if(chromos.contains(entry.chr)){
                            throw new ParserException("File is not sorted; chromosome repeated", rowCounter);
                        }
                        chromos.add(entry.chr);
                    }
                }

                indexCounter++;

                int indexSize = 10000;
                if (indexCounter == indexSize) {
                    gData.getFileIndex().add((int) reader.getPosition());
                    indexCounter = 0;
                }

                lastEntry = entry;
            }
            return gData;

        } catch (Exception e) {
            if (nextLine != null && rowCounter != 0) {
                throw new ParserException(e.getMessage(), e, rowCounter, nextLine);
            } else {
                throw new RuntimeException(e);
            }
        } finally {
            if(reader != null) reader.close();
        }
    }

    /**
     * Parse data from the given text line to {@code GWASData} instance provided
     * @param nextLine
     * @param lineNumber
     * @return  Data container, with relevant info
     * @throws ParserException If there is an error parsing the line
     *
     */
    private GWASEntry parseLine(String nextLine, long lineNumber) {
        String[] tokens = Globals.singleTabMultiSpacePattern.split(nextLine);
        if (tokens.length > 1) {

            //String chr = ParsingUtils.convertChrString(tokens[chrCol].trim());
            String chr = genome.getChromosomeAlias(tokens[this.columns.chrCol].trim());

            int start;

            try {
                start = Integer.parseInt(tokens[this.columns.locationCol].trim());
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

                return new GWASEntry(chr, start, p, nextLine);
            }
        }
        return null;
    }

    private static class GWASEntry{

        private final String chr;
        private final int start;
        private final double p;
        private final String description;

    private GWASEntry(String chr, int start, double p, String description){
        this.chr = chr;
        this.start = start;
        this.p = p;
        this.description = description;
    }
}

    /**
     * Stores numerical indexes of relevant columns
     */
    public static class GWASColumns {
        public int locationCol = -1;
        public int chrCol = -1;
        public int pCol = -1;
        public int SNPCol = -1;
        
        public boolean hasAllFields(){
            return (this.locationCol >= 0 || this.chrCol >= 0 || this.pCol >= 0 || this.SNPCol >= 0);
        }

        /**
         * Parse a columns string. Based on tokenized columns, populate column numbers to indicate which columns hold chromosome, location, p-value and SNP identifier information.
         *
         * @param headerString
         * @return
         */
        public boolean parseHeader(String headerString) {
            headerString = headerString.trim();
            String[] headers = Globals.singleTabMultiSpacePattern.split(headerString);
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

            return this.hasAllFields();
        }
    }

//    public static void main(String[] args) throws Exception{
//        generateUnsortedGWASData();
//    }

    //Used for generating test data
    public static void generateUnsortedGWASData() throws Exception{
        Random random = new Random(12345);
        int numFeats = 100;
        String headerLine = StringUtils.join(new String[]{"Chr", "bp", "p", "snp"}, "\t");

        PrintWriter writer = new PrintWriter("random.gwas");
        writer.println(headerLine);
        for(int ii=0; ii < numFeats; ii++){
            String chr = String.format("chr%d", random.nextInt(20));
            String location = String.format("%d", random.nextInt(Integer.MAX_VALUE));
            String pVal = String.format("%2.8f", random.nextDouble()/1000.0d);
            String rsID = String.format("rs%d", random.nextInt(100000));

            String line = StringUtils.join(new String[]{chr, location, pVal, rsID}, "\t");

            writer.println(line);
        }

        writer.flush();
        writer.close();

    }



}