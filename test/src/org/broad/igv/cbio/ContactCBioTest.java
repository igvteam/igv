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

package org.broad.igv.cbio;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ParsingUtils;
import org.junit.After;
import org.junit.Assume;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.List;

import static org.junit.Assert.*;

/**
 * Test loading data from cBio.
 * <p/>
 * User: jacob
 * Date: 2012/02/02
 */
public class ContactCBioTest extends AbstractHeadlessTest {

    @BeforeClass
    public static void setUpClass() throws Exception {
        AbstractHeadlessTest.setUpClass();
        checkcBioAvailable();
    }

    @After
    public void tearDown() throws Exception {
        super.tearDown();
        GeneNetwork.BASE_URL = GeneNetwork.REAL_URL;
    }

    private static void checkcBioAvailable() throws Exception {
        boolean succeeded = false;
        try {
            InputStream is = ParsingUtils.openInputStream(GeneNetwork.BASE_URL);
            is.read();
            succeeded = true;
        } catch (Exception e) {
            e.printStackTrace();
        }

        Assume.assumeTrue(succeeded);
    }

    /**
     * Load some data from cbio.
     * Checks that we are looking at the right urls
     *
     * @throws Exception
     */
    @Test
    public void testDownloadCBIO() throws Exception {
        String[] gene_list = new String[]{"egfr", "brca1", "jun"};
        GeneNetwork anno = GeneNetwork.getFromCBIO(Arrays.asList(gene_list));
        assertNotNull(anno);
    }

    /**
     * Load some data from cbio.
     * Checks that we are looking at the right urls
     *
     * @throws Exception
     */
    @Test(expected = IOException.class)
    public void testDownloadCBIOFail() throws Exception {
        String[] gene_list = new String[]{"egfr", "brca1", "jun"};
        GeneNetwork.BASE_URL += "MAKEITFAIL";
        GeneNetwork anno = GeneNetwork.getFromCBIO(Arrays.asList(gene_list));
    }

    @Test
    public void testCaching() throws Exception {
        String[] geneArray = new String[]{"sox1", "brca1", "DIRAS3"};
        List<String> geneList = Arrays.asList(geneArray);
        GeneNetwork anno = GeneNetwork.getFromCBIO(geneList);

        assertTrue(HttpUtils.isRemoteURL(anno.getSourcePath()));

        //Check that cached file exists
        String url = GeneNetwork.getURLForGeneList(geneList);
        File cachedFile = GeneNetwork.getCachedFile(url);
        assertTrue(cachedFile.exists());

        //This one should be loaded from local file
        GeneNetwork anno2 = GeneNetwork.getFromCBIO(geneList);
        assertFalse(HttpUtils.isRemoteURL(anno2.getSourcePath()));
    }

}
