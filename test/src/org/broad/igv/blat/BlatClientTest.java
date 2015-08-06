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

package org.broad.igv.blat;


import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.util.blat.BlatClient;
import org.junit.Test;

import java.util.List;

import static junit.framework.Assert.assertTrue;

/**
 * @author jrobinso
 *         Date: 11/21/12
 *         Time: 9:23 PM
 */
public class BlatClientTest extends AbstractHeadlessTest{

    /**
     * Test a basic "blat" query.
     *
     * @throws Exception
     */
    @Test
    public void testBlat() throws Exception {

        String org = "Human";
        String db = "hg19";
        String userSeq  = "AGAAGTTGGCGCAGTGGGAGACCACGTTTTATTCAGTCCAGTTCAGGATCCCCGCTATCTCAGGGCTCTCTGGGCCAGTCCTC" +
                "CTGGGAGCCCCCACCACAACACTTCCCAGGCATGAGCCCTCAGGGGCCCACATGAGCTTCCACACACTGAGAAGTGTCCGAGAAATTGGTGGG" +
                "GCCTCTGAAGGAGGCTGTGAGCAGCCCACCTGAACTCCCAGCTCACCAGCCCAAACAGGGTGCAGGGGCTCTGGCCCTGAAGAACCTGAGTGG" +
                "AGTGGAATGGCACTGGCTGGCCACTCAGCTCAGCGGGCGACGTGCCCCTACAAGTTGGCAGAAGTGGCTGCCACTGCTGGGTTTGTGTAAGAGA" +
                "GGCTGCTGCCACCATTACCTGCAGA";

        List<String> features = BlatClient.blat(org, db, userSeq);

        // We can't really assert a specific size as the blat server could get updated.  However it should have
        // several results
        assertTrue(features.size() > 1);

    }


}
