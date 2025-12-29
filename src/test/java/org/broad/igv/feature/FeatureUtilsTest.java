/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.feature;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.util.TestUtils;
import org.broad.igv.variant.vcf.VCFVariant;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.variant.vcf.VCFCodec;
import org.junit.*;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.*;

/**
 * Please create tests!!
 *
 * @author jrobinso
 */
public class FeatureUtilsTest {

    private static List<Feature> featureList;

    public FeatureUtilsTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        String dataFile = TestUtils.DATA_DIR + "bed/featureUtilsTest.bed";
        featureList = parseTestData(dataFile);

    }

    @Test
    public void testFeatureLookupInSmallVCF() throws IOException {
        String vcfFile = TestUtils.DATA_DIR + "vcf/example4-last-gsnap-2_fixed.vcf";
        Genome genome = null; // <= Don't do chromosome conversion
        FeatureCodec codec = CodecFactory.getCodec(vcfFile, genome);
        boolean isVCF = codec.getClass().isAssignableFrom(VCFCodec.class);

        TestUtils.createIndex(vcfFile);

        AbstractFeatureReader basicReader = AbstractFeatureReader.getFeatureReader(vcfFile, codec, true);
        CloseableTribbleIterator it = basicReader.iterator();
        List<VCFVariant> features = new ArrayList<VCFVariant>();
        while (it.hasNext()) {
            VCFVariant next = (VCFVariant) it.next();
            features.add(next);
        }
        // Test the first feature in the file -- this one was failing
        VCFVariant a_6321732 = (VCFVariant) FeatureUtils.getFeatureClosest(6321732.4, features);
        assertEquals("variant closest to 6321732 must be found with position=6321732", 6321732, a_6321732.getStart());


        //Now test a couple right next to each other
        VCFVariant a_6321739 = (VCFVariant) FeatureUtils.getFeatureClosest(6321739.4, features);
        VCFVariant a_6321740 = (VCFVariant) FeatureUtils.getFeatureClosest(6321740.4, features);
        assertEquals("variant closest to 6321739.4 must be found with position=6321739", 6321739, a_6321739.getStart());
        assertEquals("variant closest to 6321740.4 must be found with position=6321740", 6321740, a_6321740.getStart());
    }

    @Test
    public void testGetFeatureClosestWithIndels(){
        Feature snp = new BasicFeature("", 10, 11);
        Feature snp2 = new BasicFeature("", 15, 16);
        Feature indel = new BasicFeature("", 5, 13);
        Feature indel2 = new BasicFeature( "", 25, 31);
        List<Feature> features = List.of(indel, indel2, snp, snp2); // this order is important to make sure tiebreakers are woorking
        assertEquals(FeatureUtils.getFeatureClosest(1, features), indel);
        assertEquals(FeatureUtils.getFeatureClosest(4, features), indel);
        assertEquals(FeatureUtils.getFeatureClosest(5, features), indel);
        assertEquals(FeatureUtils.getFeatureClosest(9, features), indel);
        assertEquals("since indel and snp overlap at this point it should pick the snp", FeatureUtils.getFeatureClosest( 10, features), snp);
        assertEquals(FeatureUtils.getFeatureClosest(11, features), indel);
        assertEquals(FeatureUtils.getFeatureClosest(12, features), indel);
        assertEquals(FeatureUtils.getFeatureClosest(13, features), indel);
        assertEquals(FeatureUtils.getFeatureClosest(14, features), snp2);
        assertEquals(FeatureUtils.getFeatureClosest(15, features), snp2);
        assertEquals(FeatureUtils.getFeatureClosest(16, features), snp2);
        assertEquals(FeatureUtils.getFeatureClosest(20, features), snp2);
        assertEquals(FeatureUtils.getFeatureClosest(21, features), indel2);
        assertEquals(FeatureUtils.getFeatureClosest(100, features), indel2);
    }

    /**
     * Test of getAllFeaturesAt method, of class FeatureUtils.
     * <p/>
     * Queries the test list for features at position 500.
     */
    @Test
    public void testGetAllFeaturesAt() {

        int position = 56078756;

        List<Feature> result = FeatureUtils.getAllFeaturesContaining(position, 0, featureList);
        assertEquals(21, result.size());
        for (Feature f : result) {
            assertTrue(position >= f.getStart() && position <= f.getEnd());
        }
    }

    @Test
    public void testIndexBefore() {

        // Actual feature starts at expectedValue + 1:
        // <, 55086709, 55086713, 55955148, 56182373, >, >
        double[] positions = new double[]{1, 55086709, 55086712, 55955146, 56182372.8, 56182373, Integer.MAX_VALUE};
        int[] expectedValue = new int[]{-1, -1, 5, 29, 74, 75, 75};
        for (int i = 0; i < positions.length; i++) {
            int idx = FeatureUtils.getIndexBefore(positions[i], featureList);
            assertEquals(expectedValue[i], idx);
        }
    }

    @Test
    public void testGetFeatureAfter() {

        //chr7	55086709	55236328
        Feature f = FeatureUtils.getFeatureStartsAfter(0, featureList);
        assertEquals(55086709, f.getStart());

        //chr7	56078756	56119137
        f = FeatureUtils.getFeatureStartsAfter(56078756 - 1, featureList);
        assertEquals(56078756, f.getStart());


        //chr7	55955148	56009918
        //chr7	56019569	56024193
        f = FeatureUtils.getFeatureStartsAfter((55955148 + 56009918) / 2, featureList);
        assertEquals(56019569, f.getStart());

        // last feature
        // chr7	56182373	56184110
        f = FeatureUtils.getFeatureStartsAfter(56182373, featureList);
        assertNull(f);
    }

    @Test
    public void testGetFeatureCenterAfter() {

        //chr7	56078756	56119137
        Feature f = FeatureUtils.getFeatureCenteredAfter(56078756 - 1, featureList);
        assertEquals(56078756, f.getStart());

        //chr7	55086709	55236328
        f = FeatureUtils.getFeatureCenteredAfter(0, featureList);
        assertEquals(55086709, f.getStart());

        // last feature
        // chr7	56182373	56184110
        f = FeatureUtils.getFeatureCenteredAfter(56184110, featureList);
        assertNull(f);
    }

    @Test
    public void testGetFeatureCenterBefore() {

        // chr7 55955148	56009918
        Feature f = FeatureUtils.getFeatureCenteredBefore(56009918 + 1, featureList);
        assertEquals(55955148, f.getStart());


        f = FeatureUtils.getFeatureCenteredBefore(0, featureList);
        assertNull(f);

        // last feature
        // chr7	56182373	56184110
        f = FeatureUtils.getFeatureCenteredBefore(56184110 + 1, featureList);
        assertEquals(56182373, f.getStart());
    }

    static class TestFeature implements LocusScore {

        private int start;
        private int end;

        public TestFeature(int start, int end) {
            this.start = start;
            this.end = end;
        }


        public float getScore() {
            return 1.0f;
        }

        public LocusScore copy() {
            throw new UnsupportedOperationException("Not supported yet.");
        }

        public String getValueString(double position, int mouseX, WindowFunction windowFunction) {
            throw new UnsupportedOperationException("Not supported yet.");
        }

        public String getChr() {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        public String getContig() {
            return null;
        }

        /**
         * @return the start
         */
        public int getStart() {
            return start;
        }

        /**
         * @param start the start to set
         */
        public void setStart(int start) {
            this.start = start;
        }

        /**
         * @return the end
         */
        public int getEnd() {
            return end;
        }

        /**
         * @param end the end to set
         */
        public void setEnd(int end) {
            this.end = end;
        }
    }

    static List<Feature> parseTestData(String file) throws IOException {

        List<Feature> featureList = new ArrayList<>();
        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;
        while ((line = br.readLine()) != null) {
            String[] tokens = line.split("\t");
            featureList.add(new TestFeature(Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2])));
        }
        return featureList;

    }
}