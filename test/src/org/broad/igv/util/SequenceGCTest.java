/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
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

        String outfile = TestUtils.DATA_DIR + "/out/ecoli.wig";
        //TODO Remove the Asserts, do a check at a position instead.
        //todo This test is not very meaningful
        SequenceGC sequence = new SequenceGC(5, 1);
        sequence.ProcessPath(TestUtils.DATA_DIR + "/fasta/ecoli_out.padded.fasta.txt", outfile);

        File wigFile1 = new File(outfile);
        assertTrue(wigFile1.exists());
        assertTrue(wigFile1.getTotalSpace() > 1);
        wigFile1.deleteOnExit();
    }
}
