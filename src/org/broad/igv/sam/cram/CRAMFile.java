/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2016 UC San Diego
 * Author: Jim Robinson
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

package org.broad.igv.sam.cram;


import com.google.common.io.LittleEndianDataInputStream;
import htsjdk.samtools.cram.io.CramArray;
import htsjdk.samtools.cram.io.CramInt;
import htsjdk.samtools.cram.io.ExternalCompression;
import htsjdk.samtools.cram.io.ITF8;
import htsjdk.samtools.seekablestream.SeekableStream;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;

public class CRAMFile {

    String path;
    int major;
    int minor;
    String fileId;

    final static int RAW = 0;
    final static int GZIP = 1;
    final static int BZIP2  =2;
    final static int LZMA = 3;
    final static int RANS = 4;

    public CRAMFile(String path, int majorFormatNumber, int minorFormatNumber, String fileId) {
        this.path = path;
        this.major = majorFormatNumber;
        this.minor = minorFormatNumber;
        this.fileId = fileId;
    }


    public void readContainerHeader(long position) throws IOException {

        SeekableStream ss = IGVSeekableStreamFactory.getInstance().getStreamFor(path);
        ss.seek(position);

        BufferedInputStream bis = new BufferedInputStream(ss);

        int length = CramInt.int32(bis);

        int refSeqId = ITF8.readUnsignedITF8(bis);
        int startPos = ITF8.readUnsignedITF8(bis);
        int alignmentSpan = ITF8.readUnsignedITF8(bis);
        int nRecords = ITF8.readUnsignedITF8(bis);
        int recordCounter = ITF8.readUnsignedITF8(bis);
        int bases = ITF8.readUnsignedITF8(bis);
        int nBlocks = ITF8.readUnsignedITF8(bis);
        int[] landmarks = CramArray.array(bis);
        if (major >= 3) {
            int checksum = CramInt.int32(bis);
        }
        readBlocks(bis, nBlocks);
    }

    public void readBlocks(InputStream is, int nBlocks) throws IOException {

        LittleEndianDataInputStream lis = new LittleEndianDataInputStream(is);
        for (int i = 0; i < nBlocks; i++) {

            int compressionMethod = lis.read();
            int contentType = lis.read();
            int contentId = ITF8.readUnsignedITF8(lis);
            int size = ITF8.readUnsignedITF8(lis);
            int rawSize = ITF8.readUnsignedITF8(lis);

            byte[] blockData = new byte[size];
            lis.readFully(blockData);

            blockData = uncompress(blockData, compressionMethod);


            String tmp = new String(blockData);

            if (major >= 3) {
                int checksum = CramInt.int32(lis);
            }
        }
    }

    private byte[] uncompress(byte[] compressedContent, int method) throws IOException {

        switch (method) {
            case RAW:
                return compressedContent;
            case GZIP:
                return ExternalCompression.gunzip(compressedContent);
            case BZIP2:
                return ExternalCompression.unbzip2(compressedContent);
            case LZMA:
                return ExternalCompression.unxz(compressedContent);
            case RANS:
                return ExternalCompression.unrans(compressedContent);
            default:
                throw new RuntimeException("Unknown block compression method: " + method
                );
        }
    }
}
