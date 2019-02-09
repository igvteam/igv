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

package org.broad.igv.feature.genome.fasta;

import htsjdk.samtools.seekablestream.SeekableStream;
import org.apache.log4j.Logger;
import org.broad.igv.feature.genome.Sequence;
import org.broad.igv.util.FileUtils;
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

    private final ArrayList<String> chromoNamesList;

    public FastaIndexedSequence(String path) throws IOException {
        this(path, null);
    }

    public FastaIndexedSequence(String path, String indexPath) throws IOException {

        this.path = path;

        if (indexPath == null) indexPath = path + ".fai";

        index = new FastaIndex(indexPath);
        chromoNamesList = new ArrayList<>(index.getSequenceNames());
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
     * @param useCache
     * @return
     */

    public byte[] getSequence(String chr, int qstart, int qend, boolean useCache) {

        FastaIndex.FastaSequenceIndexEntry idxEntry = index.getIndexEntry(chr);

        if (idxEntry == null) {
            log.info("No fasta sequence entry for: " + chr);
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
            long endByte = position + endLine * bytesPerLine + offset1;

            if (startByte >= endByte) {
                return null;
            }

            // Read all the bytes in the range.  This will include endline characters
            byte[] allBytes = readBytes(startByte, endByte);

            // Create the array for the sequence -- this will be "allBytes" without the endline characters.
            ByteArrayOutputStream bos = new ByteArrayOutputStream(end - start);

            int srcPos = 0;

            // Copy first line
            final int allBytesLength = allBytes.length;
            if (offset > 0) {
                int nBases = Math.min(end - start, basesPerLine - offset);
                bos.write(allBytes, srcPos, nBases);
                srcPos += (nBases + nEndBytes);
            }

            while (srcPos < allBytesLength) {
                int nBases = Math.min(basesPerLine, allBytesLength - srcPos);
                bos.write(allBytes, srcPos, nBases);
                srcPos += (nBases + nEndBytes);
            }

            return bos.toByteArray();

        } catch (IOException e) {
            log.error("Error loading sequence " + chr + ":" + qstart + "-" + qend, e);
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
    protected byte[] readBytes(long posStart, long posEnd) throws IOException {

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

    @Override
    public boolean isRemote() {
        return FileUtils.isRemote(path);
    }

    @Override
    public boolean isFasta() {
        return true;
    }


}
