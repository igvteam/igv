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

import htsjdk.samtools.seekablestream.ISeekableStreamFactory;
import htsjdk.samtools.seekablestream.SeekableFileStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import org.broad.igv.util.HttpUtils;

import java.io.File;
import java.io.IOException;
import java.net.URL;

/**
 * @author Jim Robinson
 */
public class IGVSeekableStreamFactory implements ISeekableStreamFactory {

    private static IGVSeekableStreamFactory instance;
    static{
        instance = new IGVSeekableStreamFactory();
    }

    private IGVSeekableStreamFactory(){}

    public static IGVSeekableStreamFactory getInstance(){
        return instance;
    }

    public SeekableStream getStreamFor(URL url) throws IOException{
        return getStreamFor(url.toExternalForm());
    }

    public SeekableStream getStreamFor(String path) throws IOException {

        if (path.endsWith(".list")) {
            return new SeekableSplitStream(path);

        } else {
            SeekableStream is = null;
            if (path.toLowerCase().startsWith("http:") || path.toLowerCase().startsWith("https:")) {
                final URL url = new URL(path);
                boolean useByteRange = HttpUtils.getInstance().useByteRange(url);
                if (useByteRange) {
                    is = new IGVSeekableHTTPStream(url);
                } else {
                    is = new SeekableServiceStream(url);
                }
            } else if (path.toLowerCase().startsWith("ftp:")) {
                final URL url = new URL(path);
                is = new IGVSeekableFTPStream(url);
            } else {
                is = new SeekableFileStream(new File(path));
            }
            return is;
        }
    }

    public SeekableStream getBufferedStream(SeekableStream stream){
        return getBufferedStream(stream, IGVSeekableBufferedStream.DEFAULT_BUFFER_SIZE);
    }

    public SeekableStream getBufferedStream(SeekableStream stream, int bufferSize){
        return new IGVSeekableBufferedStream(stream, bufferSize);
    }

}
