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

package org.broad.igv.util.stream;


import org.broad.igv.util.IGVHttpUtils;
import org.broad.tribble.util.SeekableStream;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

/**
 * @author jrobinso
 * @date Jul 28, 2010
 */

public class SeekableSplitStream extends SeekableStream {

    long position = 0;
    long length = 0;
    List<PartDescriptor> descriptors;
    int currentStreamIndex = 0;

    public SeekableSplitStream(String path) throws IOException {
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
            br = org.broad.tribble.util.ParsingUtils.openBufferedReader(path);
            String nextLine;
            while ((nextLine = br.readLine()) != null) {
                String[] tokens = nextLine.split(" ");
                if (tokens.length == 2) {
                    String p = tokens[0];

                    // Require the files are in the same directory as the list file
                    String listFileName = null;
                    if(IGVHttpUtils.isURL(path)) {
                       URL url = new URL(path);
                       listFileName = (new File(url.getPath())).getName();
                    }
                    else {
                        listFileName = (new File(path)).getName();
                    }
                    SeekableStream stream = IGVSeekableStreamFactory.getStreamFor(path.replace(listFileName, p));

                    long length = Long.parseLong(tokens[1]);
                    descriptors.add(new SeekableSplitStream.PartDescriptor(length, stream));
                } else {
                    // TODO -- throw exception, or warning?
                }
            }
        }
        finally {
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
}
