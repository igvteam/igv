package org.igv.feature;

import htsjdk.tribble.Feature;
import org.igv.AbstractHeadlessTest;
import org.igv.feature.tribble.GFFCodec;
import org.igv.feature.gff.GFFFeatureSource;
import org.igv.track.TribbleFeatureSource;
import org.igv.util.ResourceLocator;
import org.igv.util.TestUtils;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import static junit.framework.Assert.*;


/**
 * @author Jim Robinson
 * @date 5/10/12
 */
public class GFFTest  extends AbstractHeadlessTest {


    private List<Feature> getFeatures(String filePath, GFFCodec.Version version) throws Exception {
        return getFeatures(filePath, "chr1", version);
    }

    private List<Feature> getFeatures(String filePath, String chr, GFFCodec.Version version) throws Exception {
        GFFFeatureSource src = new GFFFeatureSource(TribbleFeatureSource.getFeatureSource(new ResourceLocator(filePath), genome), version);
        List<Feature> features = new ArrayList<Feature>();
        Iterator<Feature> iter = src.getFeatures(chr, 0, Integer.MAX_VALUE);
        while (iter.hasNext()) features.add(iter.next());
        return features;
    }

    /**
     * Test multiple lines with the same ID, representing an alignment.
     *
     * @throws Exception
     */
    @Test
    public void testMultiLineFeature() throws Exception {
        String file = TestUtils.DATA_DIR + "gff/multi_line_feature.gff3";
        List<Feature> features = getFeatures(file, "chr1", GFFCodec.Version.GFF3);
        assertEquals(1, features.size());
        BasicFeature f = (BasicFeature) features.get(0);
        assertEquals(5, f.getExonCount());
        features = getFeatures(file, "chr2", GFFCodec.Version.GFF3);
        assertEquals(1, features.size());
        f = (BasicFeature) features.get(0);
        assertEquals(5, f.getExonCount());
    }

    /**
     * This test verifies that the attributes from column 9 are retained for a CDS feature that does not have a parent.
     * This bugfix was released with 2.1.12.
     *
     * @throws Exception
     */
    @Test
    public void testLoadSingleCDS() throws Exception {

        String gtfFile = TestUtils.DATA_DIR + "gtf/single_cds.gtf";
        List<Feature> features = getFeatures(gtfFile, GFFCodec.Version.GTF);

        assertEquals(1, features.size());
        BasicFeature bf = (BasicFeature) features.get(0);

        Exon exon = bf.getExons().get(0);
        assertNotNull(exon.getAttributes());
    }

    /**
     * Test that we parse the various spellings of "color" correctly
     * v2.2.0/2.2.1 had a bug where an attribute would be checked for
     * but a different one accessed
     *
     * @throws Exception
     */
    @Test
    public void testGFFColorSpellings() throws Exception {
        String gffFile = TestUtils.DATA_DIR + "gff/color_sps.gff";
        List<Feature> features = getFeatures(gffFile, GFFCodec.Version.GFF2);

        assertEquals("Error parsing certain features", 6, features.size());

        for (Feature f : features) {
            BasicFeature bf = (BasicFeature) f;
            boolean haveColor = bf.getColor() != null;
            if (bf.hasExons()) {
                for (Exon ex : bf.getExons()) {
                    haveColor |= ex.getColor() != null;
                }
            }
            assertTrue(haveColor);
        }
    }

    @Test
    public void testGFF2_parentIds() throws Exception {
        String path = TestUtils.DATA_DIR + "gff/parentIds.gff";
        List<Feature> features = getFeatures(path, GFFCodec.Version.GFF2);

        assertEquals(1, features.size());
        BasicFeature bf = (BasicFeature) features.get(0);

        int numExons = 17;
        assertEquals(numExons, bf.getExons().size());
        Exon lastExon = bf.getExons().get(numExons - 1);

        assertTrue(lastExon.isNonCoding());
        assertEquals(3807030 - 1, lastExon.getStart());
    }
}
