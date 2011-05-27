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
import org.broad.igv.ui.IGV;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

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
        IGV.getInstance().getGenomeManager().loadGenomeByID("hg18");
        Globals.setHeadless(true);
        dumpSequence();
    }


    /**
     * Test cached vs uncached sequence reads.
     */

    public static void dumpSequence() throws IOException {

        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("kcnn3_trf_1.fa")));

        String genome = "hg18";
        String chr = "chr1";
        int repeatStart = 153108941;
        int repeatEnd = 153108977;

        int preEnd = repeatStart;
        int preStart = preEnd - 1000;

        int postStart = repeatEnd;
        int postEnd = postStart + 1000;

        StringBuffer seqStringBuffer = new StringBuffer();
        String preSeq = new String(SequenceManager.readSequence(genome, chr, preStart, preEnd));
        String repeatString = new String(SequenceManager.readSequence(genome, chr, repeatStart, repeatEnd));
        String postSeq = new String(SequenceManager.readSequence(genome, chr, postStart, postEnd));

        seqStringBuffer.append(preSeq);
        seqStringBuffer.append(repeatString);
        seqStringBuffer.append(postSeq);
        pw.println(">kcnn3_trf_1_m" + 0);
        printFasta(pw, seqStringBuffer.toString());

        for (int i = 1; i < 5; i++) {
            pw.println(">kcnn3_trf_1_m" + i);
            seqStringBuffer = new StringBuffer();
            String modRepeat = repeatString.substring(i * 3);
            seqStringBuffer.append(preSeq);
            seqStringBuffer.append(modRepeat);
            seqStringBuffer.append(postSeq);

            String seqString = seqStringBuffer.toString();
            printFasta(pw, seqString);
        }

        String modRepeat = repeatString;
        String trf = repeatString.substring(0, 3);
        for (int i = 1; i < 10; i++) {
            pw.println(">kcnn3_trf_1_p" + i);
            seqStringBuffer = new StringBuffer();
            modRepeat = modRepeat + trf;
            seqStringBuffer.append(preSeq);
            seqStringBuffer.append(modRepeat);
            seqStringBuffer.append(postSeq);

            String seqString = seqStringBuffer.toString();
            printFasta(pw, seqString);
        }


        pw.close();
    }

    private static void printFasta(PrintWriter pw, String seqString) {
        int len = seqString.length();
        int lineLength = 80;
        for (int start = 0; start < len; start += lineLength) {
            int end = Math.min(len, start + lineLength);
            pw.println(seqString.substring(start, end));
        }
    }

}
