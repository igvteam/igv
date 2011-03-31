/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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


import org.apache.log4j.Logger;
import org.broad.tribble.util.SeekableStream;
import java.io.IOException;
import java.net.URL;

/**
 * Unfortunately the seekable stream classes exist for both Tribble and Picard, and we need both.  This class
 * is for use with Tribble and delegates all the work to a helper.
 *
 * @author jrobinso
 * @date Oct 27, 2010
 */
public class SeekableFTPStream extends SeekableStream {

    static Logger log = Logger.getLogger(SeekableFTPStream.class);

    SeekableFTPStreamHelper helper;
    public SeekableFTPStream(URL url) throws IOException {
        helper = new SeekableFTPStreamHelper(url);

    }


    public void seek(long position) {
        helper.seek(position);
    }

    public long position() {
        return helper.position();
    }

    @Override
    public boolean eof() throws IOException {
        return helper.eof();
    }

    @Override
    public long length() {
        return helper.length();
    }


    @Override
    public long skip(long n) throws IOException {
        return helper.skip(n);
    }


    @Override
    public int read(byte[] buffer, int offset, int len) throws IOException {

        return helper.read(buffer, offset, len);
    }


    public void close() throws IOException {
        helper.close();
    }

    public int read() throws IOException {
        return helper.read();
    }

}
