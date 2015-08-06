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
