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

import org.junit.BeforeClass;
import org.junit.Test;

import java.io.IOException;
import java.util.Map;

import static org.junit.Assert.assertEquals;

/**
 * @author Jim Robinson
 * @date 3/26/12
 */
public class FastaSequenceParserTest {

    static InMemorySequence fastaSequence;

    @BeforeClass
    public static void setup() throws IOException {
        String path = "http://www.broadinstitute.org/igvdata/test/fasta/ci2_test.fa";
        Map<String, byte[]> sequenceMap = FastaSequenceParser.parseFasta(path);
        fastaSequence = new InMemorySequence(sequenceMap);
    }

    @Test
    public void testReadSequence() throws Exception {


        String chr = "chr02q";
        int start = 3531385;
        int end = 3531425;

        String expectedSequence = "TAATTTTTACGTCTTATTTAAACACATATAATGAATAGGT";

        byte[] seq = fastaSequence.getSequence(chr, start, end);
        String seqString = new String(seq);

        assertEquals(expectedSequence, seqString);
    }

    @Test
    public void testReadEnd() throws Exception {

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


}


