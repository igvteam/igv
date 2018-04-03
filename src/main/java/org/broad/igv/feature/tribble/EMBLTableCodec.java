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

package org.broad.igv.feature.tribble;

import htsjdk.samtools.util.LocationAware;
import org.apache.log4j.Logger;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.Exon;
import org.broad.igv.feature.Strand;
import htsjdk.tribble.AbstractFeatureCodec;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.readers.*;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

/**
 * @author jrobinso
 *         Date: 11/15/13
 *         Time: 1:08 PM
 */
public class EMBLTableCodec extends AbstractFeatureCodec<BasicFeature, EMBLTableCodec.EmblTableIterator> {

    private static Logger log = Logger.getLogger(EMBLTableCodec.class);

    public EMBLTableCodec() {
        super(BasicFeature.class);
    }

    @Override
    public BasicFeature decode(EmblTableIterator s) throws IOException {

        EmblRecord emblRecord = s.next();
        if(emblRecord == null) return null;

        BasicFeature feature = new BasicFeature(emblRecord.getChromosome(), emblRecord.getStart(),
                emblRecord.getEnd());
        feature.setType(emblRecord.getType());
        feature.setIdentifier(emblRecord.getIdentifier());
        feature.setName(emblRecord.getIdentifier());
        feature.setStrand(emblRecord.getStrand());
        feature.setDescription(emblRecord.getDescription());
        if (emblRecord.getAlias() != null) {
            feature.setName(emblRecord.getAlias());
        }

        // If this is a "gene part" add the exons
        for (Exon exon : emblRecord.getExons()) {
            feature.addExon(exon);
        }

        return feature;
    }

    @Override
    public FeatureCodecHeader readHeader(EmblTableIterator o) throws IOException {
        return null;
    }

    /**
     * Generates a reader appropriate for use by this codec from the generic input stream.  Implementers should
     * assume the stream is buffered.
     */
    @Override
    public EmblTableIterator makeSourceFromStream(InputStream inputStream) {
        return new EmblTableIterator(inputStream);
    }

    /**
     * Generates a {@link LocationAware} reader of type EmblTableIterator.  Like {@link #makeSourceFromStream(java.io.InputStream)}, except
     * the {@link LocationAware} compatibility is required for creating indexes.
     * <p/>
     * Implementers of this method must return a type that is both {@link LocationAware} as well as EmblTableIterator.  Note that this
     * requirement cannot be enforced via the method signature due to limitations in Java's generic typing system.  Instead, consumers
     * should cast the call result into a EmblTableIterator when applicable.
     */
    @Override
    public LocationAware makeIndexableSourceFromStream(InputStream inputStream) {
        return null;    // <= indexing will fail until this is implemented.
    }

    @Override
    public boolean isDone(EmblTableIterator o) {
        return !o.hasNext();
    }

    @Override
    public void close(EmblTableIterator o) {
        o.close();
    }

    @Override
    public boolean canDecode(String path) {
        return false;
    }

    static class EmblRecord {

        private static Logger log = Logger.getLogger(EmblRecord.class);

        boolean isNegative;
        private String type;
        private String chromosome;
        private String identifier;
        private String alias;
        private String description;
        private int start = Integer.MAX_VALUE;
        private int end;
        List<Exon> exons;


        EmblRecord(String type, String chromosome, String lociString, boolean isNegative) {
            this.isNegative = isNegative;
            this.type = type;
            this.chromosome = chromosome;
            createExons(lociString, isNegative);
        }

        /**
         * Method description
         *
         * @return
         */
        public int getStart() {
            return start;
        }

        /**
         * Method description
         *
         * @return
         */
        public int getEnd() {
            return end;
        }

        /**
         * Method description
         *
         * @return
         */
        public boolean isGenePart() {
            return type.equals("CDS") || type.equals("3'UTR") || type.equals("5'UTR");
        }

        /**
         * Method description
         *
         * @return
         */
        public Strand getStrand() {
            return isNegative ? Strand.NEGATIVE : Strand.POSITIVE;
        }

        /**
         * Method description
         *
         * @return
         */
        public String getType() {
            return type;
        }

        /**
         * Method description
         *
         * @return
         */
        public String getIdentifier() {
            return identifier;
        }

        /**
         * Method description
         *
         * @param identifier
         */
        public void setIdentifier(String identifier) {
            this.identifier = identifier;
        }

        /**
         * Method description
         *
         * @return
         */
        public String getAlias() {
            return alias;
        }

        /**
         * Method description
         *
         * @param alias
         */
        public void setAlias(String alias) {
            this.alias = alias;
        }

        /**
         * Method description
         *
         * @return
         */
        public List<Exon> getExons() {
            return exons;
        }

        /**
         * Method description
         *
         * @param nextLine
         */
        public void append(String nextLine) {
            String attrString = nextLine.substring(21);
            if (attrString.startsWith("/gene=")) {
                String[] kv = attrString.split("=");
                String geneName = kv[1].replace("\"", "");
                if (geneName.startsWith("SP")) {

                    // Some genes have multiple identifiers.  Only use the first one
                    if (getIdentifier() == null) {
                        setIdentifier(geneName);
                    }
                } else {
                    setAlias(geneName);
                }
            } else if (attrString.startsWith("/systematic_id=")) {
                String[] kv = attrString.split("=");
                String id = kv[1].replace("\"", "");
                setIdentifier(id);
                setAlias(id);

            } else {
                appendToDescription(nextLine.substring(22).trim());
            }
        }

        /**
         * Method description
         *
         * @param note
         */
        public void appendToDescription(String note) {
            if (description == null) {
                description = note;
            } else {
                description += "<br>" + note;
            }
        }

        /**
         * Method description
         *
         * @return
         */
        public String getDescription() {
            return description;
        }

        /**
         * Create a list of Exon objects from the Embl join string.  Apparently exons in embl
         * format are represented by a single CDS record.
         *
         * @param joinString
         * @param isNegative
         */
        void createExons(String joinString, boolean isNegative) {
            String[] lociArray = joinString.split(",");
            exons = new ArrayList(lociArray.length);
            for (String loci : lociArray) {
                try {
                    String[] tmp = loci.split("\\.\\.");
                    int exonStart = Integer.parseInt(tmp[0]) - 1;    // - (isNegative ? 0 : 1);
                    int exonEnd = exonStart + 1;
                    if (tmp.length > 1) {
                        exonEnd = Integer.parseInt(tmp[1]);
                    }

                    Strand strand = isNegative ? Strand.NEGATIVE : Strand.POSITIVE;
                    Exon r = new Exon(chromosome, exonStart, exonEnd, strand);
                    start = Math.min(start, exonStart);
                    end = Math.max(end, exonEnd);
                    exons.add(r);


                } catch (NumberFormatException e) {
                    log.error("Error parsing exon number; " + joinString, e);
                }
            }
        }

        /**
         * Method description
         *
         * @return
         */
        public String getChromosome() {
            return chromosome;
        }
    }

    static class EmblTableIterator  {

        String chromosome;
        BufferedReader reader;
        PositionalBufferedStream is;
        EmblRecord currentRecord = null;

        EmblTableIterator(InputStream stream) {

            if(stream instanceof  PositionalBufferedStream) {
                is = (PositionalBufferedStream) stream;
            }
            else {
                is = new PositionalBufferedStream(stream);
            }

            reader = new BufferedReader(new InputStreamReader(stream));
            // Advance to first record
            next();
        }

        boolean hasNext() {
            return currentRecord != null;
        }

        EmblRecord next() {

            String nextLine = null;
            try {
                while ((nextLine = reader.readLine()) != null) {

                    if (nextLine.startsWith("ID"))    // Chromosome change
                    {
                        String chr = getFirstWord(nextLine.substring(2));
                        chromosome = chr.replace("chromosome", "chr").replace("_", "");
                    } else if (nextLine.startsWith("FT")) {
                        String featureKey = nextLine.substring(5, 19).trim();
                        if (featureKey.length() == 0) {
                            if (currentRecord != null) {
                                currentRecord.append(nextLine);
                            }
                        } else {

                            // New feature started.
                            EmblRecord returnValue = currentRecord;

                            String temp = nextLine.substring(21);
                            boolean isNegative = temp.contains("complement");
                            String lociString = parseJoinString(temp, reader).replace("<",
                                    "").replace(">", "").trim();
                            currentRecord = new EmblRecord(featureKey.trim(), chromosome, lociString, isNegative);

                            return returnValue;
                        }
                    } else {
                        // Skip line
                    }
                }

                return currentRecord;
            } catch (IOException ex) {
                log.error("Error parsing EMBL file", ex);
                return null;
            }
        }

        private String getFirstWord(String string) {
            String trimmedString = string.trim();
            char[] chars = trimmedString.toCharArray();
            int whitespaceIndex = 0;
            for (whitespaceIndex = 0; whitespaceIndex < chars.length; whitespaceIndex++) {
                if (Character.isSpaceChar(chars[whitespaceIndex])) {
                    break;
                }
            }

            return trimmedString.substring(0, whitespaceIndex).trim();

        }


        /**
         * FT   CDS             join(complement(5000933..5001976),
         * FT                   complement(5000325..5000891),complement(5000024..5000272))
         * FT                   /product="GTPase activating protein (predicted)"
         * FT                   /gene="SPAC1952.17c"
         * FT                   /gene="SPAC890.01c"
         *
         * @param joinString
         * @param reader
         * @return
         * @throws IOException
         */
        public String parseJoinString(String joinString, BufferedReader reader)
                throws IOException {

            if (joinString.startsWith("join") || joinString.startsWith("complement")) {
                int leftParenCount = countChar(joinString, '(');
                int rightParenCount = countChar(joinString, ')');
                while (leftParenCount != rightParenCount) {
                    joinString += reader.readLine().replace("FT", "").trim();
                    leftParenCount = countChar(joinString, '(');
                    rightParenCount = countChar(joinString, ')');
                }

                // join and complement functions irrelevant
                joinString = joinString.replace("join", "");
                joinString = joinString.replace("complement", "");
                joinString = joinString.replace("(", "");
                joinString = joinString.replace(")", "");
                joinString = joinString.replace('<', ' ');
                return joinString;

            } else {
                return joinString;
            }

        }


        /**
         * This must exist in the jdk ?
         *
         * @param string
         * @return
         */
        static int countChar(String string, char c) {
            int cnt = 0;
            for (int i = 0; i < string.length(); i++) {
                if (c == string.charAt(i)) {
                    cnt++;
                }
            }
            return cnt;

        }

        public void close() {
            if (reader != null) try {
                reader.close();
            } catch (IOException e) {
                log.error("Error closing EMBL reader", e);
            }
        }

    }

}