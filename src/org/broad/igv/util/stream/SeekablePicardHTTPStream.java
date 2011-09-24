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

import org.broad.tribble.util.SeekableHTTPStream;

import java.io.IOException;
import java.net.URL;

/**
 * Wrapper for the Picard interface for this class.  Its unfortunate this is required, but Picard won't include
 * Tribble and vice-versa.
 *
 * @author jrobinso
 * @date Jul 6, 2011
 */
public class SeekablePicardHTTPStream extends net.sf.samtools.util.SeekableStream {

    SeekableHTTPStream tribbleStream;
    String source;

    public SeekablePicardHTTPStream(URL url) {
        tribbleStream = new SeekableHTTPStream(new IGVUrlHelper(url));
        this.source = url.toExternalForm();
    }

    @Override
    public long length() {
       return tribbleStream.length();
    }

    @Override
    public void seek(long l) throws IOException {
        tribbleStream.seek(l);
    }

    @Override
    public int read() throws IOException {
        return tribbleStream.read();
    }

    @Override
    public int read(byte[] bytes, int offset, int length) throws IOException {
        return tribbleStream.read(bytes,  offset,  length);
    }

    @Override
    public void close() throws IOException {
        tribbleStream.close();
    }

    @Override
    public boolean eof() throws IOException {
        return tribbleStream.eof();
    }

    @Override
    public String getSource() {
        return source;
    }
}
