package org.broad.igv.feature.genome.fasta;

import org.broad.igv.feature.genome.SequenceWrapper;
import org.broad.igv.util.TestUtils;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.IOException;

import static org.junit.Assert.assertEquals;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: 8/7/11
 * Time: 9:21 PM
 * To change this template use File | Settings | File Templates.
 */
public class FastaIndexedSequenceTest {


    static FastaIndexedSequence fastaSequence;

    @BeforeClass
    public static void setup() throws IOException {
        String path = "https://s3.amazonaws.com/igv.org.test/data/ci2_test.fa";
        fastaSequence = new FastaIndexedSequence(path);
    }


    @Test
    public void testReadEnd() throws Exception {

        String path = "https://s3.amazonaws.com/igv.org.test/data/ci2_test.fa";

        String chr = "chr02q";
        int chrLen = 8059593;
        int start = chrLen - 10;
        int end = chrLen + 10;
        byte[] bytes = fastaSequence.getSequence(chr, start, end);
        assertEquals(10, bytes.length);

        byte[] expectedSequence = "TTTTTCCCAG".getBytes();

        for (int i = 0; i < 10; i++) {
            assertEquals(expectedSequence[i], bytes[i]);
        }
    }

    @Test
    public void testPaddedReference() throws Exception {

        String fasta = TestUtils.DATA_DIR + "fasta/ecoli_out.padded.fasta";
        FastaUtils.createIndexFile(fasta, fasta + ".fai");
        String expectedSequence = "atcaccattaccac******AAcggtgcgggctgacgcgtacaggaaacacagaaaaaag";
        String chr = "NC_000913_bb";
        int start = 240;
        int end = 300;
        FastaIndexedSequence sequence = new FastaIndexedSequence(fasta);

        byte[] bytes = sequence.getSequence(chr, start, end);

        assertEquals(expectedSequence, new String(bytes));

    }


    // TODO -- add some assertions, what are we testing?
    @Test
    public void testPaddedReference2() throws Exception {

        String fasta = TestUtils.DATA_DIR + "fasta/ecoli_out.padded.fasta";
        String chr = "NC_000913_bb";
        int start = 0;
        int end = 5081;
        FastaIndexedSequence sequence = new FastaIndexedSequence(fasta);


        byte[] bytes = sequence.getSequence(chr, start, end);

        for (int i = 60; i < 100; i++) {

        }

        FastaIndexedSequence s = new FastaIndexedSequence(fasta);
        SequenceWrapper sequenceHelper = new SequenceWrapper(s);
        bytes = sequenceHelper.getSequence(chr, start, end);
        for (int i = 60; i < 100; i++) {
            // ?????? Not sure what to test here
        }

    }
}
