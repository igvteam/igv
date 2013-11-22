package org.broad.igv.util.ftp;

import static junit.framework.Assert.*;

import org.junit.Test;

import java.net.URL;


/**
 * @author Jim Robinson
 * @since 10/4/11
 */
public class FTPUtilsTest {

    @Test
    public void testResourceAvailable() throws Exception {
        URL goodUrl = new URL("ftp://ftp.broadinstitute.org/pub/igv/TEST/test.txt");
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
