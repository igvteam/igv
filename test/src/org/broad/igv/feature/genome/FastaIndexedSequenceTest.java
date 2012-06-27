/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.feature.genome;

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
        String path = "http://www.broadinstitute.org/igvdata/test/fasta/ci2_test.fa";
        fastaSequence = new FastaIndexedSequence(path);
    }

    /**
     * Compare a direct read of sequence from a file vs a read from and indexed fasta for the same interval.
     *
     * @throws Exception
     */
    @Test
    public void testReadSequence() throws Exception {

        String chr02qSeqPath = "http://www.broadinstitute.org/igvdata/test/fasta/";

        String chr = "chr02q";
        int start = 3531385;
        int end = 3531425;


        // TAATTTTTACGTCTTATTTAAACACATATAATGAATAGGT;
        Sequence igvSequence = new IGVSequence(chr02qSeqPath);
        byte[] expectedBytes = igvSequence.getSequence(chr, start, end);
        String expectedSequence = new String(expectedBytes);

        byte[] bytes = fastaSequence.getSequence(chr, start, end);
        String seq = new String(bytes);
        assertEquals(expectedSequence, seq);
    }

    @Test
    public void testReadEnd() throws Exception {

        String path = "http://www.broadinstitute.org/igvdata/test/fasta/ci2_test.fa";

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
