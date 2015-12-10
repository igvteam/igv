/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2016 University of California San Diego
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
import com.mongodb.io.ByteBufferInputStream;
import htsjdk.samtools.cram.io.CramArray;
import htsjdk.samtools.cram.io.CramInt;
import htsjdk.samtools.cram.io.ITF8;
import htsjdk.samtools.seekablestream.SeekableStream;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.ByteArrayInputStream;
import java.io.IOException;

/**
 * Created by jrobinson on 12/9/15.
 */
public class CRAMFileReader {


    public static CRAMFile openFile(String path) throws IOException {

        SeekableStream stream = IGVSeekableStreamFactory.getInstance().getStreamFor(path);

        LittleEndianDataInputStream lis = new LittleEndianDataInputStream(stream);

        byte[] mn = new byte[4];
        for (int i = 0; i < 4; i++) {
            mn[i] = lis.readByte();
        }
        String magicNumber = new String(mn);

        int majorFormatNumber = lis.readUnsignedByte();
        int minorFormatNumber = lis.readUnsignedByte();

        byte[] fileId = new byte[20];

        for (int i = 0; i < 20; i++) {
            fileId[i] = lis.readByte();
        }

        String fileName = new String(fileId);

        return new CRAMFile(path, majorFormatNumber, minorFormatNumber, fileName);
    }




    public static void main(String[] args) throws IOException {

        CRAMFile cramFile = openFile(args[0]);

        cramFile.readContainerHeader(26);

    }


}
