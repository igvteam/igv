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

package org.broad.igv.util;

import junit.framework.TestCase;
import net.sf.samtools.seekablestream.SeekableHTTPStream;
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
