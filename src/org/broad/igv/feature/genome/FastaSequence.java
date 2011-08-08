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

package org.broad.igv.feature.genome;

import org.broad.igv.util.ParsingUtils;
import org.broad.tribble.util.SeekableStream;
import org.broad.tribble.util.SeekableStreamFactory;

import java.io.IOException;

/**
 * Class represents an indexed fasta file
 * Author: jrobinso
 * Date: 8/7/11
 */
public class FastaSequence {

    FastaSequenceIndex index;
    String path;
    long contentLength;

    public FastaSequence(String path) throws IOException {

        this.path = path;

        contentLength = ParsingUtils.getContentLength(path);

        // TODO -- check for existence path & index
        String indexPath = path + ".fai";
        index = new FastaSequenceIndex(indexPath);


    }


    /**
     *
     * Example:  5 bases per line, 6 bytes per line
     *
     * Bases    0 1 2 3 4 * | 5 6 7 8  9 * | 10 11 12 13 14 *  etc
     * Offset   0 1 2 3 4     0 1 2 3  4      0  1  2  3  4
     * Bytes    0 1 2 3 4 5 | 6 7 8 9 10   | 11 12 13 14 15 16
     *
     * query 9 - 13
     * start line = 1
     * base0      = 1*5 = 5
     * offset     = (9 - 5) = 4
     * start byte = (1*6) + 3 = 10
     * end   line = 2
     *
     *
     * @param chr
     * @param start
     * @param end
     * @return
     */

    public byte[] readSequence(String chr, int start, int end) {

        FastaSequenceIndex.FastaSequenceIndexEntry idxEntry = index.getIndexEntry(chr);
        if (idxEntry == null) {
            return null;
        }

        try {

            final int bytesPerLine = idxEntry.getBytesPerLine();
            final int basesPerLine = idxEntry.getBasesPerLine();
            int nEndBytes = bytesPerLine - basesPerLine;

            int startLine = start / basesPerLine;
            int endLine = end / basesPerLine;

            int base0 = startLine * basesPerLine;   // Base at beginning of start line

            int offset = start - base0;
            final long position = idxEntry.getPosition();
            long startByte = position + startLine * bytesPerLine + offset;

            int base1 = endLine * basesPerLine;
            int offset1 = end - base1;
            long endByte = Math.min(contentLength, position + endLine * bytesPerLine + offset1);

            byte[] allBytes = readBytes(startByte, endByte);

            byte[] seqBytes = new byte[end - start];

            int srcPos = 0;
            int desPos = 0;
            // Copy first line
            if(offset > 0) {
                int nBases = Math.min(end - start, basesPerLine - offset);

                System.arraycopy(allBytes, srcPos, seqBytes, desPos, nBases);
                srcPos += (nBases + nEndBytes);
                desPos += nBases;
            }

            int nBases = Math.min(basesPerLine, allBytes.length - srcPos);
            while(srcPos < allBytes.length) {
                System.arraycopy(allBytes, srcPos, seqBytes, desPos, nBases);
                srcPos += (nBases + nEndBytes);
                desPos += nBases;
            }

            return seqBytes;

        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.

            return null;
        }


    }

    private byte[] readBytes(long posStart, long posEnd) throws IOException {

        SeekableStream ss = null;
        try {
            ss = SeekableStreamFactory.getStreamFor(path);
            int nBytes = (int)(posEnd - posStart);
            byte[] bytes = new byte[nBytes];
            ss.seek(posStart);
            ss.readFully(bytes);
            return bytes;
        } finally {
            if (ss != null) {
                ss.close();
            }
        }
    }

}
