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

import net.sf.samtools.seekablestream.ISeekableStreamFactory;
import net.sf.samtools.seekablestream.SeekableFileStream;
import net.sf.samtools.seekablestream.SeekableStream;
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
