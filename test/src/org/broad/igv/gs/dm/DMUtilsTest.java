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

package org.broad.igv.gs.dm;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.PreferenceManager;
import org.broad.igv.gs.GSTestAuthenticator;
import org.broad.igv.gs.GSUtils;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.*;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.Authenticator;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static junit.framework.Assert.*;

/**
 * @author Jim Robinson
 * @date 1/16/12
 */

public class DMUtilsTest extends AbstractHeadlessTest{

    @Rule
    public TestRule testTimeout = new Timeout((int) 30e4);


    private static final String IGV_TEST_DIR = "/Home/igvtest/";
    private URL defaultURL;
    private URL personaldirectoryURL;
    private static URL fileURL;

    private static String delDirName = "testdir_deleteme";
    private static String fullPath = IGV_TEST_DIR + delDirName;

    @Before
    public void setUp() throws Exception{
        super.setUp();
        GSUtils.logout();
        //This is pretty dumb. The reason is that initializing HttpUtils sets the authenticator,
        //and we need to overwrite it in initAuth. It's not actually important which method we call,
        //as long as HttpUtils is initialized so it doesn't get initialized later
        HttpUtils.getInstance().resetAuthenticator();

        Authenticator.setDefault(new GSTestAuthenticator());

        try {
            String defaultURLStr = PreferenceManager.getInstance().get(PreferenceManager.GENOME_SPACE_DM_SERVER);
            defaultURL = new URL(defaultURLStr);
            personaldirectoryURL = new URL(defaultURLStr + DMUtils.PERSONAL_DIRECTORY);
            fileURL = new URL(defaultURL + "file");
            System.out.println("Genome space URL: " + defaultURL);
        } catch (MalformedURLException e) {
            e.printStackTrace();
            Assume.assumeTrue(false);
        }

        //We test creating this directory later
//        try{
//            DMUtils.deleteFileOrDirectory(fileURL + fullPath);
//        }catch(FileNotFoundException e){
//            //totally fine, in fact expected
//        }
    }

    @After
    public void tearDown() {
        GSUtils.logout();
        HttpUtils.getInstance().resetAuthenticator();
    }

    @Test
    public void noOpTest() {
        assertTrue(true);
    }

    @Ignore
    public void testGetDirectoryListing() throws Exception {
        final String testFileName = "Broad.080528.subtypes.seg.gz";
        boolean found = dirContainsFile(personaldirectoryURL, testFileName);
        assertTrue("Test file not found: " + testFileName, found);

    }

    protected boolean dirContainsFile(URL dirURL, String testFileName) throws Exception{
        GSDirectoryListing dirListing = DMUtils.getDirectoryListing(dirURL);

        assertNotNull("Directory listing", dirListing);

        List<GSFileMetadata> gsFiles = dirListing.getContents();

        //Search for known test file
        boolean found = false;
        for (GSFileMetadata fileMetadata : gsFiles) {
            if (fileMetadata.getName().equals(testFileName)) {
                found = true;
                String path = IGV_TEST_DIR + testFileName;
                assertEquals("Test file path not expected", path, fileMetadata.getPath());
            }
        }
        return found;
    }

    /**
     * Upload a file, check it got uploaded, delete it, check it was deleted
     * Not really ideal, but since we need to check that the file doesn't exist before
     * uploading and delete it afterwards anyway, figured we might as well combine these.
     * @throws Exception
     */
    @Ignore
    public void testUploadDeleteFile() throws Exception {

        String locName = "test2.bed";
        File localFile = new File(TestUtils.DATA_DIR + "bed", locName);
        String remPath = IGV_TEST_DIR + locName;

        // Delete, in case the file is there from a previous test run
        try {
            DMUtils.deleteFileOrDirectory(fileURL + remPath);
        } catch (IOException e) {
            // Ignore -- expected
        }

        assertFileStatus(locName, false);

        DMUtils.uploadFile(localFile, remPath);

        assertFileStatus(locName, true);

        DMUtils.deleteFileOrDirectory(fileURL + remPath);

        assertFileStatus(locName, false);
    }

    //@Test
    public void testCreateDeleteDirectory() throws Exception {

        assertFileStatus(delDirName, false);

        DMUtils.createDirectory(fileURL + fullPath);

        assertFileStatus(delDirName, true);

        DMUtils.deleteFileOrDirectory(fileURL + fullPath);

        assertFileStatus(delDirName, false);
    }

    protected void assertFileStatus(String objName, boolean expExists) throws Exception{
        boolean found = dirContainsFile(personaldirectoryURL, objName);
        if(expExists){
            assertEquals("Object not found: " + objName, expExists, found);
        }else{
            assertEquals("Object exists, but it shouldn't: " + objName, expExists, found);
        }
    }

    /**
     * Test that we can actually load files. More of an functional test than unit test.
     * We also make sure to use token authentication
     * @throws Exception
     */
    @Ignore
    public void testLoadFiles() throws Exception{

        GSDirectoryListing dirListing = DMUtils.getDirectoryListing(personaldirectoryURL);

        TrackLoader loader = new TrackLoader();
        Map<String, Exception> exceptions = new HashMap<String, Exception>();

        for(GSFileMetadata md: dirListing.getContents()){
            String mdurl = md.getUrl();
            if(!md.isDirectory() && (mdurl.endsWith(".bed") || mdurl.endsWith(".bam"))){
                System.out.println("Loading file " + mdurl);
                try{
                    List<Track> tracks = loader.load(new ResourceLocator(mdurl), genome);
                    assertNotNull(tracks);
                    assertNotSame(0, tracks.size());
                }catch(Exception e){
                    exceptions.put(mdurl, e);
                }
            }
        }

        for(Map.Entry<String, Exception> entry: exceptions.entrySet()){
            System.out.println("Exception loading Path: " + entry.getKey());
            System.out.println("StackTrace: ");
            for (StackTraceElement el : entry.getValue().getStackTrace()) {
                System.out.println(el);
            }
        }
        assertEquals(0, exceptions.size());
    }


}
