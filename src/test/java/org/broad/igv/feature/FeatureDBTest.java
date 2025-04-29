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

import htsjdk.tribble.NamedFeature;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.ui.action.SearchCommand;
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
        List<NamedFeature> features = genome.getFeatureDB().getFeaturesList(CHECK_STR, 3);
        assertEquals(3, features.size());

        features = genome.getFeatureDB().getFeaturesList(CHECK_STR, LARGE);
        assertTrue(features.size() < LARGE);
        int expected = 50;
        assertEquals(expected, features.size());
    }

    @Test
    public void testFeatureList() throws Exception {
        List<NamedFeature> features = genome.getFeatureDB().getFeaturesList(CHECK_STR, LARGE);
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
