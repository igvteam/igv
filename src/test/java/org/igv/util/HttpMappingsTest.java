package org.igv.util;

import org.junit.Test;

import java.net.MalformedURLException;

import static org.junit.Assert.*;

public class HttpMappingsTest {

    @Test
    public void mapURL() throws MalformedURLException {

        String url;
        String mappedUrl;

        // Cloudfront URL
        //https://s3.us-east-1.amazonaws.com/igv.broadinstitute.org/data/hg19/1kg/ALL.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.list
        url = "https://dn7ywbm9isq8j.cloudfront.net/data/hg19/1kg/ALL.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.list";
        mappedUrl = HttpMappings.mapURL(url);
        assertEquals("https://raw.githubusercontent.com/igvteam/igv-data/refs/heads/main/data/hg19/ALL.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.list", mappedUrl);

        // Regional URL
        url = "https://s3.us-east-1.amazonaws.com/igv.org.genomes/mm9/refGene.sorted.txt.gz";
        mappedUrl = HttpMappings.mapURL(url);
        assertEquals("https://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/refGene.txt.gz", mappedUrl);

        // Global URL
        url = "https://s3.amazonaws.com/igv.org.genomes/mm9/refGene.sorted.txt.gz";
        mappedUrl = HttpMappings.mapURL(url);
        assertEquals("https://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/refGene.txt.gz", mappedUrl);

    }
}