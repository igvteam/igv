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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.util;

import org.junit.Test;

import java.io.File;

import static org.junit.Assert.assertTrue;

public class SequenceGCTest {

    @Test
    public void testSequenceGC() throws Exception {

        String outfile = TestUtils.DATA_DIR + "out/ecoli.wig";
        //TODO Remove the Asserts, do a check at a position instead.
        //todo This test is not very meaningful
        SequenceGC sequence = new SequenceGC(5, 1);
        sequence.ProcessPath(TestUtils.DATA_DIR + "fasta/ecoli_out.padded.fasta.txt", outfile);

        File wigFile1 = new File(outfile);
        assertTrue(wigFile1.exists());
        assertTrue(wigFile1.getTotalSpace() > 1);
        wigFile1.deleteOnExit();
    }
}
