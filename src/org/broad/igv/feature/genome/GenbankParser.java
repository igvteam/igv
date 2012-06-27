package org.broad.igv.feature.genome;

import org.broad.igv.Globals;
import org.broad.igv.feature.*;
import org.broad.igv.util.ParsingUtils;
import org.broad.tribble.Feature;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Sequence defined by a Genbank (.gbk) file.  These files contain a single sequence/chromosome/contig.
 */
public class GenbankParser {

    private String accession;
    private byte[] sequence;
    private List<Feature> features;
    private String locusName;


    /**
     * @param path
     */
    GenbankParser(String path) throws IOException {

        BufferedReader reader = null;
        try {
            reader = ParsingUtils.openBufferedReader(path);
            readLocus(reader);
            readAccession(reader);
            readFeatures(reader);
            readOriginSequence(reader);
        } finally {
            if (reader != null) reader.close();
        }
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
     * @param reader
     * @throws IOException
     */
    private void readAccession(BufferedReader reader) throws  IOException {

        String line = null;
        do  {
            line = reader.readLine();
        }
        while(!line.startsWith("ACCESSION"));

        if(line == null) {
            // TODO - throw exception, missing accession
        }

        String[] tokens = Globals.whitespacePattern.split(line);
        accession = tokens[1].trim();

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

        // Skip to "FEATURES" section
        String nextLine;
        do {
            nextLine = reader.readLine();
        }
        while (!nextLine.startsWith("FEATURES"));

        //Process features until "ORIGIN"
        features = new ArrayList<Feature>();
        BasicFeature f = null;
        String currentLocQualifier = null;
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
                f.setChr(accession);
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
                                f.setAttribute(keyName, tokens[1]);
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
            int exonStart = Integer.parseInt(tmp[0]) - 1;    // - (isNegative ? 0 : 1);
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
}
