package org.igv.feature;

import htsjdk.tribble.NamedFeature;
import org.igv.AbstractHeadlessTest;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.util.List;
import java.util.Map;

import static junit.framework.Assert.*;

/**
 * User: jacob
 * Date: 2011/12/15
 */
public class FeatureDBTest extends AbstractHeadlessTest {

    public static final int LARGE = 500;

    private static final String CHECK_STR = "ABC";
    private static boolean reload = false;

    public static void main(String[] args) {


    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        AbstractHeadlessTest.setUpClass();
        reload = false;
    }

    @Before
    public void setUpTest() throws Exception {

        try {
            if (genome.getFeatureDB().size() == 0) {
                reload = true;
            }
        } catch (NullPointerException e) {
            reload = true;
        }
        if (reload) {
            setUpClass();
        }
    }

    @Test
    public void testFeaturesMap() throws Exception {
        Map<String, List<NamedFeature>> fMap = genome.getFeatureDB().getFeaturesMap(CHECK_STR);

        for (String k : fMap.keySet()) {

            assertTrue(k.startsWith(CHECK_STR));
        }

    }

    @Test
    public void testFeatureListSize() throws Exception {
        List<NamedFeature> features = genome.getFeatureDB().getFeaturesStartingWith(CHECK_STR, 3);
        assertEquals(3, features.size());

        features = genome.getFeatureDB().getFeaturesStartingWith(CHECK_STR, LARGE);
        assertTrue(features.size() < LARGE);
        int expected = 50;
        assertEquals(expected, features.size());
    }

    @Test
    public void testFeatureList() throws Exception {
        List<NamedFeature> features = genome.getFeatureDB().getFeaturesStartingWith(CHECK_STR, LARGE);
        for (NamedFeature f : features) {
            assertTrue(f.getName().startsWith(CHECK_STR));
            assertNotNull(genome.getFeatureDB().getFeature(f.getName()));
        }

    }

    @Test
    public void testMultiRetrieve() throws Exception {
        String checkstr = "EGFLAM";
        Map<String, List<NamedFeature>> fMap = genome.getFeatureDB().getFeaturesMap(checkstr);
        List<NamedFeature> data = fMap.get(checkstr);
        assertEquals(4, data.size());
    }

    @Test
    public void testMultipleEntries() throws Exception {
        String checkstr = "EG";
        Map<String, List<NamedFeature>> fMap = genome.getFeatureDB().getFeaturesMap(checkstr);
        for (String k : fMap.keySet()) {
            List<NamedFeature> data = fMap.get(k);
            //System.out.println("key " + k + " has " + data.size());
            for (int ii = 0; ii < data.size() - 1; ii++) {
                NamedFeature feat1 = data.get(ii);
                NamedFeature feat2 = data.get(ii + 1);
                int len1 = feat1.getEnd() - feat1.getStart();
                int len2 = feat2.getEnd() - feat2.getStart();
                assertTrue("Data for key " + k + " not sorted", len1 >= len2);
            }
        }
    }
}
