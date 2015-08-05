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

package org.broad.igv.feature.genome;

import htsjdk.samtools.seekablestream.SeekableStream;
import org.apache.log4j.Logger;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Implementation of Sequence backed by an indexed fasta file
 *
 * @author jrobinso
 * @date 8/7/11
 */
public class FastaIndexedSequence implements Sequence {

    static Logger log = Logger.getLogger(FastaIndexedSequence.class);

    final FastaIndex index;
    final String path;
    final long contentLength;

    private final ArrayList<String> chromoNamesList;

    public FastaIndexedSequence(String path) throws IOException {

        this.path = path;
        contentLength = ParsingUtils.getContentLength(path);

        String indexPath = path + ".fai";

// The check below is not useful in the files have been copied or moved, which is always the case for our hosted
// genomes.   It causes lots of spurious warnings
//        if(ParsingUtils.getLastModified(path) > ParsingUtils.getLastModified(indexPath)){
//            log.warn("Index file for " + path + " is older than the file it indexes");
//        }

        index = new FastaIndex(indexPath);
        chromoNamesList = new ArrayList<String>(index.getSequenceNames());
    }


    /**
     * Return the sequence for the query interval as a byte array.  Coordinates are "ucsc" style (0 based)
     * <p/>
     * Example:  5 bases per line, 6 bytes per line
     * <p/>
     * Bases    0 1 2 3 4 * | 5 6 7 8  9 * | 10 11 12 13 14 *  etc
     * Offset   0 1 2 3 4     0 1 2 3  4      0  1  2  3  4
     * Bytes    0 1 2 3 4 5 | 6 7 8 9 10   | 11 12 13 14 15 16
     * <p/>
     * query 9 - 13
     * start line = 1
     * base0      = 1*5 = 5
     * offset     = (9 - 5) = 4
     * start byte = (1*6) + 3 = 10
     * end   line = 2
     *
     * @param chr
     * @param qstart
     * @param qend
     * @return
     */

    public byte[] getSequence(String chr, int qstart, int qend) {

        FastaIndex.FastaSequenceIndexEntry idxEntry = index.getIndexEntry(chr);
        if (idxEntry == null) {
            return null;
        }

        try {

            final int start = Math.max(0, qstart);    // qstart should never be < 0
            final int end = Math.min((int) idxEntry.getSize(), qend);

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

            if (startByte >= endByte) {
                return null;
            }

            // Read all the bytes in the range.  This will include endline characters
            byte[] allBytes = readBytes(startByte, endByte);

            // Create the array for the sequence -- this will be "allBytes" without the endline characters.
            ByteArrayOutputStream bos = new ByteArrayOutputStream(end - start);

            int srcPos = 0;
            int desPos = 0;
            // Copy first line
            final int allBytesLength = allBytes.length;
            if (offset > 0) {
                int nBases = Math.min(end - start, basesPerLine - offset);
                bos.write(allBytes, srcPos, nBases);
                srcPos += (nBases + nEndBytes);
                desPos += nBases;
            }

            while (srcPos < allBytesLength) {
                int nBases = Math.min(basesPerLine, allBytesLength - srcPos);
                bos.write(allBytes, srcPos, nBases);
                srcPos += (nBases + nEndBytes);
                desPos += nBases;
            }


            return bos.toByteArray();

        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.

            return null;
        }
    }


    @Override
    public byte getBase(String chr, int position) {
        throw new RuntimeException("getBase() is not implemented for class " + FastaIndexedSequence.class.getName());
    }


    /**
     * Read the bytes between file position posStart and posEnd
     *
     * @throws IOException
     */
    private byte[] readBytes(long posStart, long posEnd) throws IOException {

        SeekableStream ss = null;
        try {
            ss = IGVSeekableStreamFactory.getInstance().getStreamFor(path);
            int nBytes = (int) (posEnd - posStart);
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

    @Override
    public List<String> getChromosomeNames() {
        return chromoNamesList;
    }

    @Override
    public int getChromosomeLength(String chrname) {
        return index.getSequenceSize(chrname);
    }
}
