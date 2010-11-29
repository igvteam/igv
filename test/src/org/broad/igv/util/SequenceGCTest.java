/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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

import static org.junit.Assert.assertEquals;
import org.junit.Test;

import java.io.File;

public class SequenceGCTest {

    @Test
    public void testSequenceGC() throws Exception{

        //TODO Remove the Asserts, do a check at a position instead.
        SequenceGC sequence = new SequenceGC(5,1);
        sequence.ProcessPath("test/data/SequenceGC/chr1.txt", "test/data/SequenceGC/chr1.wig");
        assertEquals("e145353e034285330dd9ba303b3796b3", MD5Checksum.getMD5Checksum("test/data/SequenceGC/chr1.wig"));
        File wigFile1 = new File("test/data/SequenceGC/chr1.wig");
        wigFile1.delete();
        sequence.ProcessPath("test/data/SequenceGC/chr2.txt", "test/data/SequenceGC/chr2.wig");
        assertEquals("5217101601855f1e1721b0c0fc0ea560", MD5Checksum.getMD5Checksum("test/data/SequenceGC/chr2.wig"));
        File wigFile2 = new File("test/data/SequenceGC/chr2.wig");
        wigFile2.delete();
    }
}
