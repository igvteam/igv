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

import org.broad.igv.Globals;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.StringUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.net.Authenticator;
import java.net.PasswordAuthentication;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static junit.framework.Assert.*;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: 8/8/11
 * Time: 9:16 PM
 * To change this template use File | Settings | File Templates.
 */
public class ATMUtilsTest {

    @Before
    public void setup() {
        Globals.setTesting(true);
        HttpUtils.getInstance().setAuthenticator(new GSTestAuthenticator());
    }

    @After
    public void teardown() {
        HttpUtils.getInstance().resetAuthenticator();
    }

    @Test
    public void testGetWebTools() throws Exception {
        List<WebToolDescriptor> webTools = ATMUtils.getWebTools();

        // Build a map of name -> tool
        Map<String, WebToolDescriptor> toolMap = new HashMap();
        for (WebToolDescriptor wt : webTools) {
            toolMap.put(wt.getName(), wt);
        }
        WebToolDescriptor igvDesc = toolMap.get("IGV");
        assertNotNull("IGV Descriptor", igvDesc);

    }

    @Test
    public void testGetWebTool() throws Exception {
        String toolname = "IGV";
        WebToolDescriptor igvDesc = ATMUtils.getWebTool(toolname);
        assertNotNull("IGV Descriptor", igvDesc);
        assertEquals(toolname, igvDesc.getName());

    }

    @Test
    public void testGetIGVLaunchURL() throws Exception {

        String file = "/users/igvtest/Broad.080528.subtypes.seg.gz";
        String toolname = "IGV";
        WebToolDescriptor igvDesc = ATMUtils.getWebTool(toolname);
        String url = StringUtils.decodeURL(ATMUtils.getWebtoolLaunchURL(igvDesc, file));
        assertTrue(url.startsWith(igvDesc.getBaseUrl() + "?sessionURL="));
        assertTrue(url.endsWith(file));
    }

    @Test
    public void testGetUCSCLaunchURL() throws Exception {

        List<WebToolDescriptor> webTools = ATMUtils.getWebTools();

        // Build a map of name -> tool
        Map<String, WebToolDescriptor> toolMap = new HashMap();
        for (WebToolDescriptor wt : webTools) {
            toolMap.put(wt.getName(), wt);
        }

        String toolname = ("UCSC Genome Browser");
        WebToolDescriptor ucscDesc = ATMUtils.getWebTool(toolname);

        String url = StringUtils.decodeURL(ATMUtils.getWebtoolLaunchURL(ucscDesc));
        assertNotNull(url);
    }

    static class GSTestAuthenticator extends Authenticator {

        @Override
        protected PasswordAuthentication getPasswordAuthentication() {
            return new PasswordAuthentication("igvtest", "igvtest".toCharArray());
        }
    }


}
