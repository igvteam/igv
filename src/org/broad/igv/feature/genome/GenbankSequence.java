package org.broad.igv.feature.genome;

import org.broad.igv.Globals;
import org.broad.tribble.Feature;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * Sequence defined by a Genbank (.gbk) file.  These files contain a single sequence/chromosome/contig.
 */
public class GenbankSequence implements Sequence {

    String chr;
    byte[] sequence;
    List<Feature> features;


    /**
     * Construct a sequence from the stream represented by "reader".  It is assumed that "reader" has been
     * advance to the line just after the ORIGIN keyword.
     *
     * @param reader
     */
    GenbankSequence(String chr, BufferedReader reader) throws IOException {
        this.chr = chr;
        readOriginSequence(reader);
    }

    @Override
    public byte[] readSequence(String chr, int qstart, int qend) {

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
     *
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
        do {
            nextLine = reader.readLine();

            // TODO -- first line is source (required), has total length => use to size sequence
            // TODO -- keys start at column 6,   locations and qualifiers at column 22.

            if (nextLine == null || nextLine.startsWith("ORIGIN")) {
                break;
            }

            nextLine = nextLine.trim();

            String[] tokens = Globals.whitespacePattern.split(nextLine);
            String type = tokens[0];


        } while (true);


    }


}
