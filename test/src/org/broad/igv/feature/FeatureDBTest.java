/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.feature;

import junit.framework.AssertionFailedError;
import org.broad.igv.AbstractHeadlessTest;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.util.HashMap;
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
            if (FeatureDB.size() == 0) {
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
        Map<String, List<NamedFeature>> fMap = FeatureDB.getFeaturesMap(CHECK_STR);

        for (String k : fMap.keySet()) {

            assertTrue(k.startsWith(CHECK_STR));
        }

    }

    @Test
    public void testFeatureListSize() throws Exception {
        List<NamedFeature> features = FeatureDB.getFeaturesList(CHECK_STR, 3);
        assertEquals(3, features.size());

        features = FeatureDB.getFeaturesList(CHECK_STR, LARGE);
        assertTrue(features.size() < LARGE);
        int expected = 50;
        assertEquals(expected, features.size());
    }

    @Test
    public void testFeatureList() throws Exception {
        List<NamedFeature> features = FeatureDB.getFeaturesList(CHECK_STR, LARGE);
        for (NamedFeature f : features) {
            assertTrue(f.getName().startsWith(CHECK_STR));
            assertNotNull(FeatureDB.getFeature(f.getName()));
        }

    }

    /**
     * Test thread safety by trying to read the map and clear it at the same time.
     *
     * @throws Exception
     */
    @Test
    public void testThreadSafety() throws Exception {

        final Map<Integer, AssertionFailedError> map = new HashMap<Integer, AssertionFailedError>();
        List<NamedFeature> features = FeatureDB.getFeaturesList(CHECK_STR, LARGE);
        final int expected = features.size();

        Thread read = new Thread(new Runnable() {
            public void run() {
                try {
                    List<NamedFeature> features = FeatureDB.getFeaturesList(CHECK_STR, LARGE);
                    for (NamedFeature f : features) {
                        //Check for data corruption
                        assertTrue(f.getName().startsWith(CHECK_STR));
                    }
                    assertEquals(expected, features.size());
                } catch (AssertionFailedError e) {
                    map.put(0, e);
                }
            }
        });

        Thread write = new Thread(new Runnable() {
            public void run() {
                FeatureDB.clearFeatures();
            }
        });

        read.start();
        write.start();
        read.join();

        write.join();

        features = FeatureDB.getFeaturesList(CHECK_STR, LARGE);
        assertEquals(0, features.size());

        if (map.containsKey(0)) {
            AssertionFailedError e = map.get(0);
            throw e;
        }

    }

    @Test
    public void testMultiRetrieve() throws Exception {
        String checkstr = "EGFLAM";
        Map<String, List<NamedFeature>> fMap = FeatureDB.getFeaturesMap(checkstr);
        List<NamedFeature> data = fMap.get(checkstr);
        assertEquals(4, data.size());
    }

    @Test
    public void testMultipleEntries() throws Exception {
        String checkstr = "EG";
        Map<String, List<NamedFeature>> fMap = FeatureDB.getFeaturesMap(checkstr);
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


    @Test
    public void testMutationSearch() throws Exception {

        String name = "EGFR";
        // EGFR starts with proteins MRPSG
        String[] symbols = new String[]{"M", "R", "P", "S", "G"};
        //All of these should be possible with a SNP from the EGFR sequence
        String[] muts = new String[]{"I", "R", "P", "P", "R"};
        Map<Integer, BasicFeature> matches;
        for (int ii = 0; ii < symbols.length; ii++) {
            matches = FeatureDB.getMutationAA(name, ii + 1, symbols[ii], muts[ii], genome);
            assertEquals(1, matches.size());
            for (int pos : matches.keySet()) {
                assertEquals(name, matches.get(pos).getName());
            }
        }

        name = "EGFLAM";
        int exp_start = 38439399;
        matches = FeatureDB.getMutationAA(name, 2, "H", "H", genome);
        assertEquals(1, matches.size());
        for (int geneloc : matches.keySet()) {
            assertEquals(exp_start, matches.get(geneloc).getStart());
        }

        String[] others = new String[]{"I", "M", "T"};
        for (String c : others) {
            matches = FeatureDB.getMutationAA(name, 2, "H", c, genome);
            assertEquals(0, matches.size());
        }
    }

    @Test
    public void testMutationSearchNegStrand() throws Exception {
        String name = "KRAS";
        int exp_start = 25249446;
        Map<Integer, BasicFeature> matches = FeatureDB.getMutationAA(name, 1, "M", "I", genome);
        assertEquals(1, matches.size());
        for (int geneloc : matches.keySet()) {
            assertEquals(exp_start, matches.get(geneloc).getStart());
        }

    }

    @Test
    public void testMutationSearchFail() throws Exception {
        String name = "EGFR";
        String[] symbols = "R,P,S,G,M".split(",");
        Map<Integer, BasicFeature> matches;
        for (int ii = 0; ii < symbols.length; ii++) {
            matches = FeatureDB.getMutationAA(name, ii + 1, symbols[ii], "M", genome);
            assertEquals(0, matches.size());
        }
    }

    @Test
    public void testMutationSearchNT() throws Exception {
        String name = "EGFR";
        String[] bps = new String[]{"A", "T", "G"};
        Map<Integer, BasicFeature> matches;
        for (int ii = 0; ii < bps.length; ii++) {
            matches = FeatureDB.getMutationNT(name, ii + 1, bps[ii], genome);
            assertEquals(1, matches.size());
        }
    }

    @Test
    public void testMutationSearchNTNegStrand() throws Exception {
        String name = "KRAS";
        String[] bps = new String[]{"A", "T", "G"};
        Map<Integer, BasicFeature> matches;
        for (int ii = 0; ii < bps.length; ii++) {
            matches = FeatureDB.getMutationNT(name, ii + 1, bps[ii], genome);
            assertEquals(1, matches.size());
        }

        //Exon 3 of KRAS, starting at amino acid 56
        int startNT = (56 - 1) * 3;
        char[] bps2 = "CTCGACACAGCAGGT".toCharArray();
        for (int ii = 0; ii < bps2.length; ii++) {
            matches = FeatureDB.getMutationNT(name, ii + startNT + 1, String.valueOf(bps2[ii]), genome);
            assertEquals(1, matches.size());
        }
    }

}
