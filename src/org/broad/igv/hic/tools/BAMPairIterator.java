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

package org.broad.igv.hic.tools;

import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.ReadMate;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;

import java.io.IOException;
import java.util.Map;

/**
 * @author Jim Robinson
 * @date 9/24/11
 */
public class BAMPairIterator implements PairIterator {

    AlignmentPair nextPair = null;
    CloseableIterator<Alignment> iterator;
    private AlignmentReader reader;
    // Map of name -> index
    private Map<String, Integer> chromosomeOrdinals;

    public BAMPairIterator(String path) throws IOException {

        this.reader = AlignmentReaderFactory.getReader(path, false);
        this.iterator = reader.iterator();
        advance();
    }

    private void advance() {

        while (iterator.hasNext()) {
            Alignment alignment = iterator.next();

            final ReadMate mate = alignment.getMate();
            if (alignment.isPaired() && alignment.isMapped() && alignment.getMappingQuality() > 0 &&
                    mate != null && mate.isMapped()) {
                // Skip "normal" insert sizes
                if ((!alignment.getChr().equals(mate.getChr())) || alignment.getInferredInsertSize() > 1000) {

                    // Each pair is represented twice in the file,  keep the record with the "leftmost" coordinate

                    if ((alignment.getChr().equals(mate.getChr()) && alignment.getStart() < mate.getStart()) ||
                            (alignment.getChr().compareTo(mate.getChr()) < 0)) {
                        final String chrom1 = alignment.getChr();
                         final String chrom2 = mate.getChr();
                        if (chromosomeOrdinals.containsKey(chrom1) && chromosomeOrdinals.containsKey(chrom2)) {
                            int chr1 = chromosomeOrdinals.get(chrom1);
                            int chr2 = chromosomeOrdinals.get(chrom2);
                           //  nextPair = new AlignmentPair(chr1, alignment.getStart(), chr2, mate.getStart());
                        }
                        return;
                    }
                }
            }

        }
        nextPair = null;

    }

    public boolean hasNext() {
        return nextPair != null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public AlignmentPair next() {
        AlignmentPair p = nextPair;
        advance();
        return p;
    }

    public void remove() {
        // Not implemented
    }

    public void close() {
        iterator.close();
        try {
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

}
