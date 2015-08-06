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

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.util.TestUtils;
import org.junit.Before;
import org.junit.Test;

import java.io.*;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.Set;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: 8/7/11
 * Time: 7:41 PM
 * To change this template use File | Settings | File Templates.
 */
public class FastaIndexTest extends AbstractHeadlessTest {

    static String indexPath = TestUtils.DATA_DIR + "fasta/ci2_test.fa.fai";
    private FastaIndex index;

    public static void main(String[] args) throws IOException {
        long numContigs = (long) 1.4e6;
        String outfilepath = TestUtils.DATA_DIR + numContigs + "contigs.fasta";
        genManyContigFasta(outfilepath, 80, numContigs);
    }

    /**
     * Generate a fasta file with a large number of contigs
     * Each will be only 1 line
     *
     * @param outfilepath
     * @param linelen
     * @param numContigs
     */
    static void genManyContigFasta(String outfilepath, int linelen, long numContigs) throws IOException {
        char[] chars = {'A', 'C', 'G', 'T'};
        long seed = 3458056478l;
        Random rand = new Random(seed);
        BufferedWriter out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outfilepath)));
        PrintWriter writer = new PrintWriter(out);
        String line;

        for (int num = 0; num < numContigs; num++) {

            writer.print(">contig");
            writer.println(num);

            line = "";
            for (int ch = 0; ch < linelen; ch++) {
                line += chars[rand.nextInt(chars.length)];
            }
            writer.println(line);
        }

        writer.flush();
        writer.close();

    }

    @Before
    public void setUp() throws IOException {
        index = new FastaIndex(indexPath);
    }

    @Test
    public void testGetContigs() throws Exception {

        List expectedContigs = Arrays.asList("chr01p", "chr02q", "chr03q");

        Set<String> contigs = index.getSequenceNames();
        assertEquals(expectedContigs.size(), contigs.size());
        for (Object exp : expectedContigs) {
            assertTrue(contigs.contains(exp));
        }

    }

    @Test
    public void testGetIndexEntry() throws Exception {
        //5803340	14565324	50	51
        FastaIndex.FastaSequenceIndexEntry entry = index.getIndexEntry("chr03q");
        assertEquals("Size", 5803340, entry.getSize());
        assertEquals("Position", 14565324, entry.getPosition());
        assertEquals("basesPerLine", 50, entry.getBasesPerLine());
        assertEquals("bytesPerLine", 51, entry.getBytesPerLine());
    }


}
