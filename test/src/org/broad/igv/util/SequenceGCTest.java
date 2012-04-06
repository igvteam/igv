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
