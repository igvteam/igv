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

package org.igv.util;

import org.junit.Test;

import java.io.IOException;

import static junit.framework.Assert.assertFalse;
import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

/**
 * @author Jim Robinson
 * @date Jul 27, 2011
 */

public class IGVHttpClientUtilsTest {

    String url = "https://s3.amazonaws.com/igv.org.test/csi_test/chr10p.bam.csi";
    int byteCount = 105063;

    @Test
    public void testGetContentLength() throws IOException {

        assertEquals(byteCount, HttpUtils.getInstance().getContentLength(HttpUtils.createURL(url)));
    }

    @Test
    public void testExists() throws IOException {

        assertTrue("Resource unexpectedly does not exist", HttpUtils.getInstance().resourceAvailable(url));

        url = "http://nosuchserver/genomes/hg18.genome";
        assertFalse("Resource unexpectedly found", HttpUtils.getInstance().resourceAvailable(url));

        url = "http://igvdata.broadinstitute.org/nosuchfile.txt";
        assertFalse(HttpUtils.getInstance().resourceAvailable(url));
    }

}
