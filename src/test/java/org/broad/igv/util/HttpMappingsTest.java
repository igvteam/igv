package org.broad.igv.util;

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
        assertEquals("https://raw.githubusercontent.com/igvteam/igv-genomes/refs/heads/main/data/hg19/ALL.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.list", mappedUrl);

        // Regional URL
        url = "https://s3.us-east-1.amazonaws.com/igv.org.genomes/hg38/rmsk/hg38_rmsk_DNA.bed.gz";
        mappedUrl = HttpMappings.mapURL(url);
        assertEquals("https://hgdownload.soe.ucsc.edu/hubs/GCF/000/001/405/GCF_000001405.40/bbi/GCF_000001405.40_GRCh38.p14.rmsk.DNA.bb", mappedUrl);

        // Global URL
        url = "https://s3.amazonaws.com/igv.org.genomes/hg38/rmsk/hg38_rmsk_DNA.bed.gz";
        mappedUrl = HttpMappings.mapURL(url);
        assertEquals("https://hgdownload.soe.ucsc.edu/hubs/GCF/000/001/405/GCF_000001405.40/bbi/GCF_000001405.40_GRCh38.p14.rmsk.DNA.bb", mappedUrl);


    }
}