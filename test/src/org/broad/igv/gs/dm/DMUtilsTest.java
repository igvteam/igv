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

package org.broad.igv.gs.dm;

import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.gs.GSUtils;
import org.broad.igv.util.HttpUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.net.Authenticator;
import java.net.PasswordAuthentication;
import java.net.URL;
import java.util.List;

import static junit.framework.Assert.*;

/**
 * @author Jim Robinson
 * @date 1/16/12
 */
public class DMUtilsTest {


    @Before
    public void setup() {
        Globals.setTesting(true);
        GSUtils.logout();
        HttpUtils.getInstance().setAuthenticator(new GSTestAuthenticator());
    }

    @After
    public void teardown() {
        //       HttpUtils.getInstance().resetAuthenticator();
    }

    @Test
    public void testListPersonalDirectory() throws Exception {

        URL defaultURL = new URL(PreferenceManager.getInstance().get(PreferenceManager.GENOME_SPACE_DM_SERVER) +
                DMUtils.PERSONAL_DIRECTORY);
        GSDirectoryListing dirListing = DMUtils.getDirectoryListing(defaultURL);
        assertNotNull("Directory listing", dirListing);

        List<GSFileMetadata> gsFiles = dirListing.getContents();

        //Search for known test file
        boolean found = false;
        final String testFileName = "Broad.080528.subtypes.seg.gz";
        for (GSFileMetadata fileMetadata : gsFiles) {
            if (fileMetadata.getName().equals(testFileName)) {
                found = true;
                String path = "/Home/igvtest/" + testFileName;
                assertEquals("Test file path", path, fileMetadata.getPath());
            }
        }
        if (!found) {
            fail("Test file not found: " + testFileName);
        }
    }

    @Test
    public void testGetDirectoryListing() throws Exception {
        // TODO -- test for arbitrary directory
    }

    @Test
    public void testUploadFile() throws Exception {
        // TODO -- implementation
    }

    @Test
    public void testCreateDirectory() throws Exception {
        // TODO -- implementation
    }


    static class GSTestAuthenticator extends Authenticator {

        @Override
        protected PasswordAuthentication getPasswordAuthentication() {
            return new PasswordAuthentication("igvtest", "igvtest".toCharArray());
        }
    }

}
