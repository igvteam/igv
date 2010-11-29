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
package org.broad.igv.maf;

//~--- non-JDK imports --------------------------------------------------------

import edu.mit.broad.prodinfo.genomicplot.ParseException;
import edu.mit.broad.prodinfo.multiplealignment.MAFAlignment;
import edu.mit.broad.prodinfo.multiplealignment.MAFIO;
import edu.mit.broad.prodinfo.multiplealignment.MultipleAlignment.AlignedSequence;
import edu.mit.broad.prodinfo.sequence.Sequence;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * @author jrobinso
 */
public class VirusMafTest {

    /**
     * Method description
     *
     * @param args
     * @throws IOException
     * @throws ParseException
     */
    public static void main(String[] args) throws Exception {
//chr21:39,511,925-39,511,951
        String in = "/Users/jrobinso/IGV/virus/poxvirus.1.USA_2003_044.maf";
        printSeqNames(in);
        load(in);
        //loadTiles(in);
    }

    private static void printSeqNames(String in) throws IOException, ParseException {
        MAFIO mafio = new MAFIO();

        MAFAlignment maf = mafio.load(in);

        System.out.println(maf.getReferenceId());

        for (AlignedSequence seq : maf.getAlignedSequences()) {
            System.out.println(seq.getId());
        }


    }

    private static void load(String in) throws ParseException, IOException {

        // chr21:34,854,786-34,854,830
        int start = 0;
        int end = 100;

        MAFIO mafio = new MAFIO();

        List<String> seqToLoad = new ArrayList();
        seqToLoad.add("USA_2003_044");
        seqToLoad.add("USA_2003_039");
        MAFAlignment maf = mafio.load(in, seqToLoad, start, end);


        maf.load(in, start, end);

        List<AlignedSequence> sequences = maf.getAlignedSequences();
        System.out.println("Ref start = " + maf.getReferenceStart());

        AlignedSequence seq0 = sequences.get(0);

        for (Integer i : seq0.getGapSizes()) {
            System.out.print(i + " ");
        }
        System.out.println();
        System.out.println("Ungapped length = " + seq0.getUngappedLength());
        System.out.println("Seq length = " + seq0.getLength());
        System.out.println("Seq size =" + seq0.getSequenceBases().length());
        System.out.println("Gap=" + seq0.getGapsSize() + "\t");
        for (int i = 0; i < (seq0.getLength()); i++) {
            System.out.print(seq0.getGapAdjustedCoordinate(i) + " ");
        }
        System.out.println();
        for (int i = 0; i < (seq0.getLength()); i++) {
            System.out.print(seq0.getSequencePosition(i) + " ");
        }
        System.out.println();

        for (AlignedSequence sequence : sequences) {
            Sequence seq = sequence.getSequence();

            System.out.println("Start = " + sequence.getStart());

            for (int i = 0; i < (seq.getSequenceBases().length()); i++) {
                System.out.print(seq0.getGapAdjustedCoordinate(i) + " ");
            }
            System.out.println();

            //System.out.print("Gap="+seq.getGapsSize() + "\t");

            System.out.print(sequence.getId() + "\t");
            System.out.println(seq.getSequenceBases());

        }
    }
}
