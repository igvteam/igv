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

package org.broad.igv.util;

import junit.framework.TestCase;
import htsjdk.samtools.seekablestream.SeekableHTTPStream;
import org.broad.igv.util.stream.IGVSeekableHTTPStream;
import org.broad.igv.util.stream.SeekableServiceStream;
import org.junit.Test;

import java.net.URL;

/**
 * @author jrobinso
 * @date Jul 28, 2010
 */
public class SeekableServiceStreamTest extends TestCase {

    /**
     * Test a file at some random position using the webservice, and compare results obtained to the standard
     * http stream.
     *
     * @throws Exception
     */
    @Test
    public void testRead() throws Exception {

        String tdfFile = "http://www.broadinstitute.org/igvdata/annotations/hg18/conservation/omega.12mer.tdf";

        HttpUtils.getInstance().updateProxySettings();

        IGVSeekableHTTPStream hs = new IGVSeekableHTTPStream(new URL(tdfFile));
        final int position = 100;
        hs.seek(position);
        final int range = 1000;
        byte[] expectedBytes = new byte[range];
        hs.read(expectedBytes, 0, expectedBytes.length);

        SeekableServiceStream sss = new SeekableServiceStream(new URL(tdfFile));
        sss.seek(position);
        byte[] bytes = new byte[range];
        sss.read(bytes, 0, bytes.length);

        for (int i = 0; i < expectedBytes.length; i++) {
            assertEquals(expectedBytes[i], bytes[i]);
        }
    }


}
