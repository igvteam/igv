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

package org.broad.igv.blat;


import org.broad.igv.AbstractHeadlessTest;
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

        List<String []> features = BlatClient.blat(org, db, userSeq);

        // We can't really assert a specific size as the blat server could get updated.  However it should have
        // several results
        assertTrue(features.size() > 1);

    }


}
