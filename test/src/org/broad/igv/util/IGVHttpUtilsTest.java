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
import org.broad.igv.Globals;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: May 26, 2010
 * Time: 1:46:38 PM
 * To change this template use File | Settings | File Templates.
 */
public class IGVHttpUtilsTest extends TestCase {

    public void testGetETag() throws Exception {

        String url = "http://www.broadinstitute.org/igvdata/test/NA12878.SLX.bam";

        String eTag = IGVHttpUtils.getETag(new URL(url));

        assertEquals("\"47c51-2aa59e7866-486ff338f6d35\"", eTag);
    }


    public void testGetContentLength() throws IOException {

        String url = "http://igvdata.broadinstitute.org/genomes/hg18.genome";
        assertEquals(3617657, IGVHttpUtils.getContentLength(new URL(url)));
    }

    public void testExists() throws IOException {
        URL url = new URL("http://igvdata.broadinstitute.org/genomes/hg18.genome");
        assertTrue(IGVHttpUtils.resourceAvailable(url));

        url = new URL("http://nosuchserver/genomes/hg18.genome");
        assertFalse(IGVHttpUtils.resourceAvailable(url));

        url = new URL("http://igvdata.broadinstitute.org/nosuchfile.txt");
        assertFalse(IGVHttpUtils.resourceAvailable(url));
    }
}
