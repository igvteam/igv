/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.sam.reader;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.CloseableIterator;
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
