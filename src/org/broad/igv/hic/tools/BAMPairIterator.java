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
import org.broad.igv.sam.reader.AlignmentQueryReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.Iterator;

/**
 * @author Jim Robinson
 * @date 9/24/11
 */
public class BAMPairIterator implements PairIterator {

    AlignmentPair nextPair = null;
    CloseableIterator<Alignment> iterator;
    private AlignmentQueryReader reader;

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
                    nextPair = new AlignmentPair(alignment.getChr(), alignment.getStart(), mate.getChr(),
                            mate.getStart());
                    return;
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
