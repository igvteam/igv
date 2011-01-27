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

package org.broad.igv.util;

import org.broad.igv.Globals;
import org.broad.igv.feature.SequenceManager;
import org.broad.igv.feature.genome.GenomeManager;

import java.io.IOException;

import static org.junit.Assert.assertEquals;

/**
 * A one-off class to generate alternative reference sequences around a repeat region.  Checked into SVN so I don't
 * loose it.
 *
 * @author jrobinso
 * @date Jan 27, 2011
 */
public class AltSequenceGenerator {


    public static void main(String[] args) throws IOException {
        GenomeManager.getInstance().findGenomeAndLoad("hg18");
        Globals.setHeadless(true);
        dumpSequence();
    }


    /**
     * Test cached vs uncached sequence reads.
     */

    public static void dumpSequence() {

        String genome = "hg18";
        String chr = "chr1";
        int repeatStart = 153108941;
        int repeatEnd = 153108977;

        int preEnd = repeatStart;
        int preStart = preEnd - 10;

        int postStart = repeatEnd;
        int postEnd = postStart + 10;

        byte[] preSeq = SequenceManager.readSequence(genome, chr, preStart, preEnd);
        String seqString = new String(preSeq);
        System.out.println(seqString);

        byte[] repeatSeq = SequenceManager.readSequence(genome, chr, repeatStart, repeatEnd);
        seqString = new String(repeatSeq);
        System.out.println(seqString);

        byte[] postSeq = SequenceManager.readSequence(genome, chr, postStart, postEnd);
        seqString = new String(postSeq);
        System.out.println(seqString);

    }

}
