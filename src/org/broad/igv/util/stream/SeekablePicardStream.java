/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.util.stream;


import net.sf.samtools.seekablestream.SeekableStream;

import java.io.BufferedInputStream;
import java.io.IOException;

/**
 * Unfortunately the seekable stream classes exist for both Tribble and Picard, and we need both.  This class
 * is for use with Picard (SAM) and delegates all the work to a helper.
 *
 * @author jrobinso
 * @date Dec 23, 2011
 * @deprecated Should no longer be necessary since tribble and picard were merged
 */
@Deprecated
public class SeekablePicardStream extends net.sf.samtools.seekablestream.SeekableStream {

    public static final int DEFAULT_BUFFER_SIZE = 1024000;

    final private int bufferSize;
    SeekableStream wrappedStream;
    String source;
    BufferedInputStream bufferedStream;
    long position;


    public SeekablePicardStream(SeekableStream tribbleStream, String source) {

        this.wrappedStream = tribbleStream;
        this.source = source;
        this.position = 0;
        this.bufferSize = DEFAULT_BUFFER_SIZE;
        bufferedStream = new BufferedInputStream(wrappedStream, bufferSize);
    }

    @Override
    public long length() {
        return wrappedStream.length();
    }

    @Override
    public long position() throws IOException {
        return this.position;
    }

    @Override
    public void seek(long position) throws IOException {
        this.position = position;
        wrappedStream.seek(position);
        bufferedStream = new BufferedInputStream(wrappedStream, bufferSize);
    }

    @Override
    public int read() throws IOException {
        int b = bufferedStream.read();
        position++;
        return b;
    }


    @Override
    public int read(byte[] buffer, int offset, int length) throws IOException {
        int nBytesRead = bufferedStream.read(buffer, offset, length);
        if (nBytesRead > 0) {
            position += nBytesRead;
        }
        return nBytesRead;
    }


    @Override
    public void close() throws IOException {
        wrappedStream.close();
    }

    @Override
    public boolean eof() throws IOException {
        return wrappedStream.eof();
    }

    @Override
    public String getSource() {
        return source;
    }
}
