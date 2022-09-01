package org.broad.igv.feature.genome.fasta;

import org.broad.igv.feature.genome.Sequence;

import java.io.IOException;

public class CompareSeqs {

    //static String fasta1="http://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg38/hg38.fa";
    static String fasta1 = "/Users/jrobinso/Downloads/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta";
    static String fasta2 = "/Users/jrobinso/Downloads/Homo_sapiens_assembly38/GRCh38_full_analysis_set_plus_decoy_hla.fa";
    //static String compressedSequencePath="https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg38/hg38.fa.gz";

    public static void main(String [] args) throws IOException {

        //Sequence fastaSequence = new FastaIndexedSequence(sequenceLocation);
        Sequence seq_1 = new FastaIndexedSequence(fasta1);
        Sequence seq_2 = new FastaIndexedSequence(fasta2); //new FastaBlockCompressedSequence(compressedSequencePath);

        int len1 = seq_1.getChromosomeLength("chr12");
        int len2 = seq_2.getChromosomeLength("chr12");

        if (len1 != len2) {
            System.err.println("Chr lengths not equal");
        }

        byte[] seq1 = seq_1.getSequence("chr12", 0, len1, false);
        byte[] seq2 = seq_2.getSequence("chr12", 0, len1, false);

        for (int i = 0; i < len1; i++) {
            byte b1 = seq1[i];
            if (b1 >= 97) b1 -= 32;

            byte b2 = seq2[i];
            if (b2 >= 97) b2 -= 32;

            if (b1 != b2) {
                System.err.println("Seq mismatch at position " + (i + 1) + "   " + Character.toString(b1) + " <> " + Character.toString(b2));
                return;
            }
        }
        System.out.println("Success");
    }
}
