/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTIES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
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

import org.junit.Test;

import java.io.File;
import java.net.URL;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * @author Jim Robinson
 * @date 1/17/12
 */
public class HttpUtilsTest {
    @Test
    public void testGetContentsAsString() throws Exception {

    }

    @Test
    public void testGetContentLength() throws Exception {

    }

    @Test
    public void testCompareResources() throws Exception {

        File localFile = new File(TestUtils.DATA_DIR + "/bed/test.bed");
        assertTrue(localFile.exists());

        // To make this test robust we must compare the local file time to the known remote one
        long remoteModifedTime = 1326815079000l;
        boolean expectedValue = localFile.lastModified() >= remoteModifedTime;

        URL remoteURL = new URL("http://www.broadinstitute.org/igvdata/test/bed/test.bed");
        boolean comp = HttpUtils.getInstance().compareResources(localFile, remoteURL);

        assertEquals(expectedValue, comp);
    }


}
