package org.broad.igv.gs.dm;

import org.broad.igv.Globals;
import org.broad.igv.gs.GSUtils;
import org.broad.igv.util.HttpUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.net.Authenticator;
import java.net.PasswordAuthentication;
import java.util.List;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertNotNull;
import static junit.framework.Assert.fail;

/**
 * @author Jim Robinson
 * @date 1/16/12
 */
public class DMUtilsTest {


    @Before
    public void setup() {
        Globals.setTesting(true);
       // HttpUtils.getInstance().setAuthenticator(new GSTestAuthenticator());
        GSUtils.clearGSToken();
    }

    @After
    public void teardown() {
 //       HttpUtils.getInstance().resetAuthenticator();
    }

    @Test
    public void testListDefaultDirectory() throws Exception {

        GSDirectoryListing dirListing = DMUtils.listDefaultDirectory();
        assertNotNull("Directory listing", dirListing);

        List<GSFileMetadata> gsFiles = dirListing.getContents();

        //Search for known test file
        boolean found = false;
        final String testFileName = "Broad.080528.subtypes.seg.gz";
        for (GSFileMetadata fileMetadata : gsFiles) {
            if (fileMetadata.getName().equals(testFileName)) {
                found = true;
                String path = "/users/igvtest/" + testFileName;
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
