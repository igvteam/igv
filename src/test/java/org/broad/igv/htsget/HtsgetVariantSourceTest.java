package org.broad.igv.htsget;

import htsjdk.tribble.Feature;
import htsjdk.variant.vcf.VCFHeader;
import org.broad.igv.feature.genome.Genome;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Iterator;

import static org.junit.Assert.*;

@Ignore   // Problems with reference server?
public class HtsgetVariantSourceTest {

    @Test
    public void testGetHeader() throws Exception {

        String url = "https://htsget.ga4gh.org/variants/giab.NA12878";
        Genome genome = null;
        HtsgetUtils.Metadata metadata =  HtsgetUtils.getMetadata(url);
        HtsgetVariantSource source = new HtsgetVariantSource(metadata, genome);
        VCFHeader header = (VCFHeader) source.getHeader();
        assertNotNull(header);
    }

    @Test
    public void testReadFeatures() throws Exception {

        String url = "https://htsget.ga4gh.org/variants/giab.NA12878";
        String chr = "8";
        int start = 128732400 - 1;
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
        assertEquals(11, featureCount);


    }

//    public static void main(String [] args) throws Exception {
//        (new HtsgetVariantSourceTest()).testReadFeatures();
//    }

}