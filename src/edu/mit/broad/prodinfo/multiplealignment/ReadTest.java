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

package edu.mit.broad.prodinfo.multiplealignment;

import edu.mit.broad.prodinfo.multiplealignment.*;
import edu.mit.broad.prodinfo.genomicplot.ParseException;
import edu.mit.broad.prodinfo.multiplealignment.MultipleAlignment.AlignedSequence;
import edu.mit.broad.prodinfo.sequence.Sequence;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


/**
*
* @author jrobinso
*/
public class ReadTest {
    
    public static void main(String [] args) throws IOException, ParseException {
        
        String in = "/Users/jrobinso/IGV/MAF/chr21.maf";
        
        //chr21:34,854,786-34,854,830
        int start = 34854773;
        int end = 34854841 + 1000;
        
        MAFIO mafio = new MAFIO();
        
        MAFAlignment maf = mafio.load(in, null, start, end) ;
        
        long t0 = System.currentTimeMillis();
        List<AlignedSequence> sequences = maf.getAlignedSequences();
        long dt = System.currentTimeMillis() - t0;
        System.out.println("Elapsed time = " + dt);
        
        for(AlignedSequence sequence : sequences) {
            Sequence seq = sequence.getSequence();
            System.out.print(sequence.getName() + "\t");
            System.out.println(seq.getSequenceBases());
        }

    }

}

