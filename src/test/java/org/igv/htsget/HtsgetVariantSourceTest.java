package org.igv.htsget;

import htsjdk.tribble.Feature;
import htsjdk.variant.vcf.VCFHeader;
import org.igv.feature.genome.Genome;
import org.junit.Test;

import java.util.Iterator;

import static org.junit.Assert.*;

public class HtsgetVariantSourceTest {

    @Test
    public void testGetHeader() throws Exception {

        String url = "https://htsget.ga4gh-demo.org/variants/spec-v4.3";
        Genome genome = null;
        HtsgetUtils.Metadata metadata =  HtsgetUtils.getMetadata(url);
        HtsgetVariantSource source = new HtsgetVariantSource(metadata, genome);
        VCFHeader header = (VCFHeader) source.getHeader();
        assertNotNull(header);
    }

    @Test
    public void testReadFeatures() throws Exception {

        String url = "https://htsget.ga4gh-demo.org/variants/spec-v4.3";
        String chr = "20";
        int start = 0;
        int end = 128770475;
        Genome genome = null;

        HtsgetUtils.Metadata metadata =  HtsgetUtils.getMetadata(url);
        HtsgetVariantSource source = new HtsgetVariantSource(metadata, genome);
        Iterator<Feature> featureIterator = source.getFeatures(chr, start, end);
        int featureCount = 0;
        while (featureIterator.hasNext()) {
            Feature f = featureIterator.next();
            featureCount++;
        }
        assertTrue( featureCount > 0);


    }

//    public static void main(String [] args) throws Exception {
//        (new HtsgetVariantSourceTest()).testReadFeatures();
//    }

}