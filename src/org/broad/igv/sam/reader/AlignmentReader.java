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

package org.broad.igv.sam.reader;

import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.sam.Alignment;

import java.io.IOException;
import java.util.List;
import java.util.Set;

/**
 * @author jrobinso
 */
public interface AlignmentReader {

    void close() throws IOException;

    /**
     * Return the list of sequence (chromosome) names as defined in the files header or meta-data section.
     */
    List<String> getSequenceNames();

    /**
     * Return the set of all platforms represented in this file.
     * May return "null"
     */
    Set<String>  getPlatforms();

    CloseableIterator<Alignment> iterator();

    /**
     * Query alignments over a given range. Be careful about start/end,
     * SAMTools uses 1-based, but IGV uses 0-based.
     * This function requires hasIndex() == true.
     *
     *
     * @param sequence
     * @param start 0-based start location
     * @param end 0-based, exclusive-end coordinate
     * @param contained
     * @return
     * @throws IOException
     */
    CloseableIterator<Alignment> query(final String sequence, final int start, final int end, final boolean contained) throws IOException;

    boolean hasIndex();

}
