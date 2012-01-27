/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.feature;

import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.util.TestUtils;
import org.broad.igv.variant.vcf.VCFVariant;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.iterators.CloseableTribbleIterator;
import org.broad.tribble.source.BasicFeatureSource;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;
import org.junit.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.*;

/**
 * Please create tests!!
 *
 * @author jrobinso
 */
public class FeatureUtilsTest {

    static List<Feature> features;

    public FeatureUtilsTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {

        features = new ArrayList();

        for (int i = 0; i < 1000; i++) {
            features.add(new TestFeature(i, i + 5));
        }

    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    @Test
    public void testFeatureLookupInSmallVCF() throws IOException {
        String vcfFile = TestUtils.DATA_DIR + "/vcf/example4-last-gsnap-2_fixed.vcf";
        FeatureCodec codec = CodecFactory.getCodec(vcfFile);
        boolean isVCF = codec.getClass().isAssignableFrom(VCFCodec.class);
        BasicFeatureSource basicReader = BasicFeatureSource.getFeatureSource(vcfFile, codec, true);
        CloseableTribbleIterator it = basicReader.iterator();
        List<VCFVariant> features = new ArrayList<VCFVariant>();
        while (it.hasNext()) {
            VCFVariant next = (VCFVariant) it.next();
            features.add(next);
        }
        //System.out.println("Now looking up closest to..");
        VCFVariant a_6321739 = (VCFVariant) FeatureUtils.getFeatureClosest(6321739d, features);
        VCFVariant a_6321740 = (VCFVariant) FeatureUtils.getFeatureClosest(6321740d, features);
        //System.out.println(a_6321739);
        //System.out.println(a_6321739);
        assertEquals("variant closest to 6321739 must be found with position=6321739", 6321739, a_6321739.getStart());
        assertEquals("variant closest to 6321740 must be found with position=6321740", 6321740, a_6321740.getStart());


    }

    /**
     * Test of divideByChromosome method, of class FeatureUtils.
     */
    //@Test
    public void testDivideByChromosome() {
        System.out.println("divideByChromosome");
        List<IGVFeature> features = null;
        Map<String, List<IGVFeature>> expResult = null;
        Map<String, List<IGVFeature>> result = FeatureUtils.divideByChromosome(features);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of segreateFeatures method, of class FeatureUtils.
     */
    //@Test
    public void testSegreateFeatures() {
        System.out.println("segreateFeatures");
        List<IGVFeature> features = null;
        double scale = 0.0;
        List<List<IGVFeature>> expResult = null;
        List<List<IGVFeature>> result = FeatureUtils.segreateFeatures(features, scale);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of sortFeatureList method, of class FeatureUtils.
     */
    //@Test
    public void testSortFeatureList() {
        System.out.println("sortFeatureList");
        List<? extends LocusScore> features = null;
        FeatureUtils.sortFeatureList(features);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getFeatureAt method, of class FeatureUtils.
     */
    //@Test
    public void testGetFeatureAt() {
        System.out.println("getFeatureAt");
        double position = 0.0;
        double minWidth = 0.0;
        List<? extends LocusScore> features = null;
        LocusScore expResult = null;
        LocusScore result = (LocusScore) FeatureUtils.getFeatureAt(position, 0, features);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getAllFeaturesAt method, of class FeatureUtils.
     * <p/>
     * Queries the test list for features at position 500.
     */
    @Test
    public void testGetAllFeaturesAt() {
        System.out.println("getFeatureAt");
        int position = 500;
        int maxLength = 100;

        List<Feature> result = FeatureUtils.getAllFeaturesAt(position, maxLength, 0, features);
        assertEquals(6, result.size());
        for (Feature f : result) {
            assertTrue(position >= f.getStart() && position <= f.getEnd());
        }


    }

    /**
     * Test of getIndexAfter method, of class FeatureUtils.
     */
    @Test
    public void testGetIndexAfter() {

        List<LocusScore> features = new ArrayList();
        for (int i = 0; i <= 10000; i += 5) {
            int start = i;
            int end = start + 10;
            features.add(new TestFeature(start, end));
        }


        int expResult = 99;
        int result = FeatureUtils.getIndexBefore(499, features);
        assertEquals(expResult, result);

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

        public String getValueString(double position, WindowFunction windowFunction) {
            throw new UnsupportedOperationException("Not supported yet.");
        }

        public String getChr() {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
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
}