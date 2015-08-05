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

package org.broad.igv.feature.genome;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.*;
import org.broad.igv.util.ParsingUtils;
import htsjdk.tribble.Feature;
import org.broad.igv.util.StringUtils;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Sequence defined by a Genbank (.gbk) file.  These files contain a single sequence/chromosome/contig.
 */
public class GenbankParser {

    private static Logger log = Logger.getLogger(GenbankParser.class);

    private String path;
    private String accession;
    private byte[] sequence;
    private List<Feature> features;
    private String locusName;
    private String[] aliases;

    private static List<String> nameFields = Arrays.asList("gene");


    /**
     * @param path
     */
    public GenbankParser(String path) throws IOException {

        this.path = path;
        readFeatures(true);
    }

    public GenbankParser() throws IOException {

    }


    public void readFeatures(boolean readSequence) throws IOException {
        BufferedReader reader = null;
        try {
            reader = ParsingUtils.openBufferedReader(path);
            readFeatures_(readSequence, reader);
        } finally {
            if (reader != null) reader.close();
        }
    }

    public void readFeatures(InputStream inputStream, boolean readSequence) throws IOException {
        BufferedReader reader = null;
        try {
            reader = new BufferedReader(new InputStreamReader(inputStream));
            readFeatures_(readSequence, reader);
        } finally {
            if (reader != null) reader.close();
        }
    }

    private void readFeatures_(boolean readSequence, BufferedReader reader) throws IOException {
        readLocus(reader);

        String line = null;
        do {
            line = reader.readLine();
            if (line.startsWith("ACCESSION")) {
                readAccession(line);
            } else if (line.startsWith("ALIASES")) {
                readAliases(line);
            }
        }
        while (line != null && !line.startsWith("FEATURES"));

        readFeatures(reader);
        if (readSequence) readOriginSequence(reader);
    }


    public byte[] getSequence(String chr, int qstart, int qend) {

        if (sequence == null) {
            return null;
        } else {
            final int start = Math.max(0, qstart);    // qstart should never be < 0
            final int end = Math.min(sequence.length, qend);
            int len = end - start;

            byte[] bytes = new byte[len];
            Arrays.fill(bytes, (byte) 0);
            int s = Math.max(start, 0);
            System.arraycopy(sequence, s, bytes, 0, len);  // This copy is a crime
            return bytes;
        }
    }


    /**
     * Return the total sequence length.  Added primarily for unit testing.
     *
     * @return
     */
    public int getSequenceLenth() {
        return sequence == null ? 0 : sequence.length;
    }


    /**
     * Read the locus line
     * LOCUS       NT_030059             105338 bp    DNA     linear   CON 28-OCT-2010
     */
    private void readLocus(BufferedReader reader) throws IOException {
        String line = reader.readLine();
        String[] tokens = Globals.whitespacePattern.split(line);
        if (!tokens[0].equalsIgnoreCase("LOCUS")) {
            // throw exception
        }
        locusName = tokens[1].trim();
    }

    /**
     * Read the acession line
     * ACCESSION   K03160
     *
     * @throws IOException
     */
    private void readAccession(String line) {

        String[] tokens = Globals.whitespacePattern.split(line);
        if (tokens.length < 2) {
            log.info("Genbank file missing ACCESSION number.");
        } else {
            accession = tokens[1].trim();
        }
    }


    /**
     * Read the sequence aliases line  -- Note: this is an IGV extension
     * ACCESSION   K03160
     *
     * @throws IOException
     */
    private void readAliases(String line) {
        String[] tokens = Globals.whitespacePattern.split(line);
        if (tokens.length < 2) {
            //log.info("Genbank file missing ACCESSION number.");
        } else {
            aliases = Globals.commaPattern.split(tokens[1]);
        }
    }


    /**
     * Read the origin section.   Example...
     * <p/>
     * ORIGIN
     * 1 gatcctccat atacaacggt atctccacct caggtttaga tctcaacaac ggaaccattg
     * 61 ccgacatgag acagttaggt atcgtcgaga gttacaagct aaaacgagca gtagtcagct
     * 121 ctgcatctga agccgctgaa gttctactaa gggtggataa catcatccgt gcaagaccaa
     *
     * @param reader
     */
    private void readOriginSequence(BufferedReader reader) throws IOException {

        String nextLine;
        ByteArrayOutputStream buffer = new ByteArrayOutputStream(100000); // TODO -- we might know length

        while ((nextLine = reader.readLine()) != null && !nextLine.startsWith("//")) {
            nextLine = nextLine.trim();
            String[] tokens = Globals.whitespacePattern.split(nextLine);
            for (int i = 1; i < tokens.length; i++) {
                buffer.write(tokens[i].getBytes());
            }
        }
        sequence = buffer.toByteArray();
    }


    /**
     * Return a string representing the chromosome/contig/sequence.  We use the accession if it is defined, otherwise
     * the first word in the LOCUS field.
     *
     * @return
     */
    public String getChr() {
        return accession == null ? locusName : accession;
    }

    /**
     * FEATURES             Location/Qualifiers
     * source          1..105338
     * /organism="Homo sapiens"
     * /mol_type="genomic DNA"
     * /db_xref="taxon:9606"
     * /chromosome="10"
     * gene            1..105338
     * /gene="PTEN"
     * /note="Derived by automated computational analysis using
     * gene prediction method: BestRefseq."
     * /db_xref="GeneID:5728"
     * /db_xref="HGNC:9588"
     * /db_xref="HPRD:03431"
     * /db_xref="MIM:601728"
     * <p/>
     * CDS             join(1033..1111,30588..30672,62076..62120,67609..67652,
     * 69576..69814,88681..88822,94416..94582,97457..97681,
     * 101850..102035)
     * /gene="PTEN"
     *
     * @param reader
     * @throws IOException
     */
    private void readFeatures(BufferedReader reader) throws IOException {


        String chr = getChr();

        //Process features until "ORIGIN"
        features = new ArrayList<Feature>();
        BasicFeature f = null;
        String currentLocQualifier = null;
        String nextLine = null;
        do {
            nextLine = reader.readLine();

            // TODO -- first line is source (required), has total length => use to size sequence
            // TODO -- keys start at column 6,   locations and qualifiers at column 22.

            if (nextLine == null || nextLine.startsWith("ORIGIN")) {
                break;
            }

            if (nextLine.charAt(5) != ' ') {
                String featureType = nextLine.substring(5, 21).trim();
                f = new BasicFeature();
                f.setChr(chr);
                f.setType(featureType);
                currentLocQualifier = nextLine.substring(21);

                if (!featureType.toLowerCase().equals("source")) {
                    features.add(f);
                }


            } else {
                String tmp = nextLine.substring(21).trim();
                if (tmp.length() > 0)
                    if (tmp.charAt(0) == '/') {
                        if (currentLocQualifier.charAt(0) == '/') {
                            String[] tokens = Globals.equalPattern.split(currentLocQualifier, 2);
                            if (tokens.length > 1) {
                                String keyName = tokens[0].length() > 1 ? tokens[0].substring(1) : "";
                                String value = StringUtils.stripQuotes(tokens[1]);
                                f.setAttribute(keyName, value);
                                if (nameFields.contains(keyName)) {
                                    f.setName(value);
                                }
                            } else {
                                // TODO -- don't know how to interpret, log?
                            }
                        } else {
                            // location string TODO -- many forms of this to support
                            // Crude test for strand
                            Strand strand = currentLocQualifier.contains("complement") ? Strand.NEGATIVE : Strand.POSITIVE;
                            f.setStrand(strand);


                            // join and complement functions irrelevant
                            String joinString = currentLocQualifier.replace("join", "");
                            joinString = joinString.replace("order", "");
                            joinString = joinString.replace("complement", "");
                            joinString = joinString.replace("(", "");
                            joinString = joinString.replace(")", "");

                            if (joinString.contains("..")) {

                                joinString = joinString.replace("<", "");
                                joinString = joinString.replace(">", "");

                                List<Exon> exons = createExons(joinString, strand);
                                FeatureUtils.sortFeatureList(exons);
                                Exon firstExon = exons.get(0);
                                f.setStart(firstExon.getStart());
                                Exon lastExon = exons.get(exons.size() - 1);
                                f.setEnd(lastExon.getEnd());
                                if (exons.size() > 1) {
                                    for (Exon exon : exons) {
                                        f.addExon(exon);
                                    }
                                }
                            } else {
                                // TODO Single locus for now,  other forms possible
                                int start = Integer.parseInt(joinString) - 1;
                                int end = start + 1;
                                f.setStart(start);
                                f.setEnd(end);
                            }

                        }
                        currentLocQualifier = tmp;
                    } else {
                        currentLocQualifier = (currentLocQualifier == null ? tmp : currentLocQualifier + tmp);
                    }
            }
        }
        while (true);
    }


    /**
     * Create a list of Exon objects from the Embl join string.  Apparently exons in embl
     * format are represented by a single CDS record.
     *
     * @param joinString
     */
    List<Exon> createExons(String joinString, Strand strand) {
        String[] lociArray = joinString.split(",");
        List<Exon> exons = new ArrayList(lociArray.length);
        boolean isNegative = joinString.contains("complement");
        for (String loci : lociArray) {
            String[] tmp = loci.split("\\.\\.");
            int exonStart = 0;    // - (isNegative ? 0 : 1);
            try {
                exonStart = Integer.parseInt(tmp[0]) - 1;
            } catch (NumberFormatException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
            int exonEnd = exonStart + 1;
            if (tmp.length > 1) {
                exonEnd = Integer.parseInt(tmp[1]);
            }

            Exon r = new Exon(accession, exonStart, exonEnd, strand);
            exons.add(r);
        }
        return exons;
    }

    public String getAccession() {
        return accession;
    }

    public byte[] getSequence() {
        return sequence;
    }

    public List<Feature> getFeatures() {
        return features;
    }

    public String getLocusName() {
        return locusName;
    }

    public String[] getAliases() {
        return aliases;
    }
}
