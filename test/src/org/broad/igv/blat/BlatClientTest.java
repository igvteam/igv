package org.broad.igv.blat;


import org.broad.igv.blat.BlatClient;
import org.junit.Test;

import java.util.List;

import static junit.framework.Assert.assertTrue;

/**
 * @author jrobinso
 *         Date: 11/21/12
 *         Time: 9:23 PM
 */
public class BlatClientTest {

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
