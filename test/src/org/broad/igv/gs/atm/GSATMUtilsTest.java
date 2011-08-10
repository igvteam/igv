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

package org.broad.igv.gs.atm;

import org.broad.igv.gs.dm.DMUtils;
import org.broad.igv.util.BrowserLauncher;
import org.junit.Test;

import java.net.URL;
import java.net.URLEncoder;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: 8/8/11
 * Time: 9:16 PM
 * To change this template use File | Settings | File Templates.
 */
public class GSATMUtilsTest {

    static final String BASE_URL = "https://atmtest.genomespace.org:8443/atm/webtools";


    @Test
    public void testGetWebTools() throws Exception {
        List<WebToolDescriptor> webTools = ATMUtils.getWebTools(new URL(BASE_URL));
        for (WebToolDescriptor wt : webTools) {
            wt.print();
        }

    }

    @Test
    public void testGetWebtoolLaunchURL() throws Exception {

        String url = ATMUtils.getWebtoolLaunchURL(BASE_URL, "Cytoscape");
        System.out.println(url);

    }


    @Test
    public void testGetSubtoolLaunchURL() throws Exception {
        URL defaultURL = new URL(DMUtils.baseUrl + "defaultdirectory");
        System.out.println(DMUtils.getDirectoryListing(defaultURL).getDirectory());
        String gsPath = URLEncoder.encode("/users/test/259.wgs.muc1.hist.txt", "UTF-8");
        String url = ATMUtils.getSubtoolLaunchURL(BASE_URL, "GenePattern", "ConvertLineEndings?input.filename=" + gsPath);
        System.out.println(url);

        BrowserLauncher.openURL(url);

    }
}
