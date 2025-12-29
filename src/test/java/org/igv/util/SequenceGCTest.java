/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.igv.util;

import org.junit.Test;

import java.io.File;

import static org.junit.Assert.assertTrue;

public class SequenceGCTest {

    @Test
    public void testSequenceGC() throws Exception {

        String outfile = TestUtils.TMP_OUTPUT_DIR + "ecoli.wig";
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
