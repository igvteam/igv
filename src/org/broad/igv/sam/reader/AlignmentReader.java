/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.sam.reader;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.sam.Alignment;

import java.io.IOException;
import java.util.List;
import java.util.Set;

/**
 * @author jrobinso
 */
public interface AlignmentReader<T extends Alignment> {

    void close() throws IOException;

    /**
     * Return the list of sequence (chromosome) names as defined in the files header or meta-data section.
     */
    List<String> getSequenceNames();

    /**
     * Return the header of the SAM file. May be null
     * @return
     */
    SAMFileHeader getFileHeader();

    /**
     * Return the set of all platforms represented in this file.
     * May return "null"
     */
    Set<String>  getPlatforms();

    CloseableIterator<T> iterator();

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
    CloseableIterator<T> query(final String sequence, final int start, final int end, final boolean contained) throws IOException;

    boolean hasIndex();

}
