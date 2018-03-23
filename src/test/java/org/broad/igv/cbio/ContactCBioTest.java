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

package org.broad.igv.cbio;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ParsingUtils;
import org.junit.*;

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
@Ignore
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
     * Load some data from cBio.
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
