/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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
