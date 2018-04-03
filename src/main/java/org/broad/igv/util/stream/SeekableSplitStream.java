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

package org.broad.igv.util.stream;


import htsjdk.samtools.seekablestream.SeekableStream;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

/**
 * Class to support a single logical file split into multiple parts.  Introduced to support Amazon cloud files that
 * might be split, but not used thus far.
 *
 * @author jrobinso
 * @date Jul 28, 2010
 */

public class SeekableSplitStream extends SeekableStream {

    long position = 0;
    long length = 0;
    List<PartDescriptor> descriptors;
    int currentStreamIndex = 0;

    private String path;

    public SeekableSplitStream(String path) throws IOException {
        this.path = path;
        parseDescriptors(path);
    }

    public void seek(long position) throws IOException {
        this.position = position;
        long end = 0;
        long start = 0;
        for (PartDescriptor desc : descriptors) {
            end += desc.getContentLength();
            if (position >= start && position < end) {
                long delta = position - start;
                desc.getStream().seek(delta);
            } else {
                desc.getStream().seek(0);
            }
            start = end;
        }
    }

    public long position() throws IOException {
        return position;
    }

    @Override
    public int read(byte[] buffer, int off, int len) throws IOException {

        int bytesRead = 0;
        long end = 0;
        long start = 0;
        for (PartDescriptor desc : descriptors) {
            end += desc.getContentLength();
            if (position >= start && position < end) {
                int delta = (int) (desc.contentLength - desc.stream.position());
                int l = Math.min(len - bytesRead, delta);
                int nBytes = desc.stream.read(buffer, off + bytesRead, l);
                if (nBytes < 0) {
                    return nBytes;
                }
                bytesRead += nBytes;
                position += nBytes;
                if (bytesRead >= len) {
                    break;
                }
            }
            start = end;
        }


        // If we get this far, and haven't read any bytes, we're at EOF
        return bytesRead > 0 ? bytesRead : -1;
    }

    @Override
    public int read() throws IOException {
        int b = descriptors.get(currentStreamIndex).getStream().read();
        position++;
        return b;
    }

    @Override
    public void close() throws IOException {
        for (PartDescriptor desc : descriptors) {
            desc.getStream().close();
        }
    }

    private void parseDescriptors(String path) throws IOException {

        BufferedReader br = null;
        descriptors = new ArrayList();
        try {
            br = new BufferedReader(new InputStreamReader(ParsingUtils.openInputStream(path)));
            String nextLine;
            while ((nextLine = br.readLine()) != null) {
                String[] tokens = nextLine.split(" ");
                if (tokens.length == 2) {
                    String p = tokens[0];

                    // Require the files are in the same directory as the list file
                    String listFileName = null;
                    if (HttpUtils.isRemoteURL(path)) {
                        URL url = new URL(path);
                        listFileName = (new File(url.getPath())).getName();
                    } else {
                        listFileName = (new File(path)).getName();
                    }
                    SeekableStream stream = IGVSeekableStreamFactory.getInstance().getStreamFor(path.replace(listFileName, p));

                    long length = Long.parseLong(tokens[1]);
                    descriptors.add(new SeekableSplitStream.PartDescriptor(length, stream));
                } else {
                    // TODO -- throw exception, or warning?
                }
            }
        } finally {
            if (br != null) {
                br.close();
            }
        }
    }


    public long length() {
        return length;
    }


    public static class PartDescriptor {
        private long contentLength;
        private SeekableStream stream;

        public PartDescriptor(long contentLength, SeekableStream stream) {
            this.contentLength = contentLength;
            this.stream = stream;
        }

        public long getContentLength() {
            return contentLength;
        }

        public SeekableStream getStream() {
            return stream;
        }
    }

    @Override
    public boolean eof() throws IOException {
        return false;  //TODO -- punting on this for the moment
    }

    @Override
    public String getSource() {
        return this.path;
    }
}
