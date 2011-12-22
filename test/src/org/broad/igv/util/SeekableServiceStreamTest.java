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

package org.broad.igv.util;

import junit.framework.TestCase;
import org.broad.igv.util.stream.IGVUrlHelper;
import org.broad.igv.util.stream.SeekableServiceStream;
import org.broad.tribble.util.SeekableHTTPStream;
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

        SeekableHTTPStream hs = new SeekableHTTPStream(new IGVUrlHelper(new URL(tdfFile)));
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
