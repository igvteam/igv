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
import org.junit.Test;

import java.io.IOException;
import java.net.URL;

import static junit.framework.Assert.assertFalse;
import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: May 26, 2010
 * Time: 1:46:38 PM
 * To change this template use File | Settings | File Templates.
 */
public class IGVHttpClientUtilsTest  {

    @Test
    public void testGetContentLength() throws IOException {

        String url = "http://igvdata.broadinstitute.org/genomes/hg18.genome";
        assertEquals(3617644, IGVHttpClientUtils.getContentLength(new URL(url)));
    }

    @Test
    public void testExists() throws IOException {
        URL url = new URL("http://igvdata.broadinstitute.org/genomes/hg18.genome");
        assertTrue(IGVHttpClientUtils.resourceAvailable(url));

        url = new URL("http://nosuchserver/genomes/hg18.genome");
        assertFalse(IGVHttpClientUtils.resourceAvailable(url));

        url = new URL("http://igvdata.broadinstitute.org/nosuchfile.txt");
        assertFalse(IGVHttpClientUtils.resourceAvailable(url));
    }
}
