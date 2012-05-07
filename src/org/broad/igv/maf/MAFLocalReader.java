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
package org.broad.igv.maf;

import edu.mit.broad.prodinfo.genomicplot.ParseException;
import edu.mit.broad.prodinfo.multiplealignment.MAFAlignment;
import edu.mit.broad.prodinfo.multiplealignment.MAFIO;
import edu.mit.broad.prodinfo.multiplealignment.MultipleAlignment.AlignedSequence;
import org.apache.log4j.Logger;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 */
public class MAFLocalReader implements MAFReader {

    private static Logger log = Logger.getLogger(MAFLocalReader.class);
    public static String mafFile;
    MAFIO mafIO;
    List<String> sequenceIds;

    public MAFLocalReader(String mafFileName) throws IOException, ParseException {
        mafFile = mafFileName;
        mafIO = new MAFIO(mafFile, false);
    }

    /**
     * Return the sequence names, as in species (not chromosomes).
     *
     * @return
     */
    public List<String> getChrNames() {
        return sequenceIds;
    }

    public MAFTile loadTile(String seq, int start, int end, List<String> species) {
        try {

            MAFAlignment maf = mafIO.load(species, start, end);

            if (maf.isEmpty()) {
                return null;
            }

            AlignedSequence refSeq = maf.getAlignedSequences().get(0);

            // TODO -- compare refSeq with seq, if not equal => error

            int[] gapAdjustedCoordinates = new int[end - start];
            for (int i = start; i < end; i++) {
                int idx = i - start;
                gapAdjustedCoordinates[idx] = refSeq.getGapAdjustedCoordinate(idx);
            }

            Map<String, String> bases = new HashMap();
            for (String seqId : maf.getAlignedSequenceIds()) {
                bases.put(seqId, maf.getAlignedSequence(seqId).getSequenceBases());
            }

            return new MAFTile(start, end, bases, gapAdjustedCoordinates, refSeq.getId());

        } catch (Exception ex) {
            ex.printStackTrace();
            return null;
        }
    }

    @Override
    protected void finalize() throws Throwable {
        super.finalize();
        mafIO.destroyFileHandle();
    }
}
