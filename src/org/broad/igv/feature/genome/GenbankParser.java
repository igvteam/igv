package org.broad.igv.feature.genome;

import org.broad.igv.Globals;
import org.broad.igv.feature.*;
import org.broad.igv.util.StringUtils;
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
public class GenbankParser  {

    String chr;
    byte[] sequence;
    List<Feature> features;


    /**
     * Construct a sequence from the stream represented by "reader".  It is assumed that "reader" has been
     * advance to the line just after the ORIGIN keyword.
     *
     * @param reader
     */
    GenbankParser(String chr, BufferedReader reader) throws IOException {
        this.chr = chr;
        readFeatures(reader);
        readOriginSequence(reader);
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
                                f.setAttribute(tokens[0], tokens[1]);
                            } else {
                                // TODO -- don't know how to interpret, log?
                            }
                        } else {
                            // location string TODO -- many forms of this to support

                            if (currentLocQualifier.contains ("..")) {
                                List<Exon> exons = parseJoinString(currentLocQualifier, chr);
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
                            }
                            else {
                                // TODO Single locus for now,  other forms possible
                                int start = Integer.parseInt(currentLocQualifier) - 1;
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
     * FT   CDS             join(complement(5000933..5001976),
     * FT                   complement(5000325..5000891),complement(5000024..5000272))
     * FT                   /product="GTPase activating protein (predicted)"
     * FT                   /gene="SPAC1952.17c"
     * FT                   /gene="SPAC890.01c"
     *
     * @param joinString
     * @param chr
     * @return
     * @throws IOException
     */
    public static List<Exon> parseJoinString(String joinString, String chr)
            throws IOException {

        if (joinString.startsWith("join") || joinString.startsWith("complement")) {
            int leftParenCount = StringUtils.countChar(joinString, '(');
            int rightParenCount = StringUtils.countChar(joinString, ')');
            while (leftParenCount != rightParenCount) {
                leftParenCount = StringUtils.countChar(joinString, '(');
                rightParenCount = StringUtils.countChar(joinString, ')');
            }

            // join and complement functions irrelevant
            joinString = joinString.replace("join", "");
            joinString = joinString.replace("complement", "");
            joinString = joinString.replace("(", "");
            joinString = joinString.replace(")", "");
            joinString = joinString.replace('<', ' ');
            return createExons(joinString, chr);

        } else {
            return createExons(joinString, chr);
        }

    }


    /**
     * Create a list of Exon objects from the Embl join string.  Apparently exons in embl
     * format are represented by a single CDS record.
     *
     * @param joinString
     * @param chromosome
     */
    static List<Exon> createExons(String joinString, String chromosome) {
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

            Strand strand = isNegative ? Strand.NEGATIVE : Strand.POSITIVE;
            Exon r = new Exon(chromosome, exonStart, exonEnd, strand);
            exons.add(r);
        }
        return exons;
    }
}
