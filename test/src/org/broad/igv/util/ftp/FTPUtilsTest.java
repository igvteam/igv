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

package org.broad.igv.util.ftp;

import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.net.URL;

import static junit.framework.Assert.*;


/**
 * @author Jim Robinson
 * @since 10/4/11
 */
public class FTPUtilsTest {

    @Test
    public void testResourceAvailable() throws Exception {
        URL goodUrl = new URL(TestUtils.AVAILABLE_FTP_URL);
        assertTrue(FTPUtils.resourceAvailable(goodUrl));

        URL nonExistentURL = new URL("ftp://ftp.broadinstitute.org/pub/igv/TEST/doesntExist");
        assertFalse(FTPUtils.resourceAvailable(nonExistentURL));

        URL nonExistentServer = new URL("ftp://noSuchServer/pub/igv/TEST/doesntExist");
        assertFalse(FTPUtils.resourceAvailable(nonExistentServer));
    }

    @Test
    public void testNCBIResourceAvailable() throws Exception {

        URL vcfURL = new URL("ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz");
        assertTrue(FTPUtils.resourceAvailable(vcfURL));

        URL tabixURL = new URL("ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz.tbi");
        assertTrue(FTPUtils.resourceAvailable(tabixURL));

    }

    @Test
    public void testGetContentLength() throws Exception {

        URL vcfURL = new URL("ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz");

        long contentLength = 1297662912;

        assertEquals(contentLength, FTPUtils.getContentLength(vcfURL));

    }
}
