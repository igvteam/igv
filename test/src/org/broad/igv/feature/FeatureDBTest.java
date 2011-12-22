/*
 * Copyright (c) 2007-2012 by The Broad Institute of MIT and Harvard.All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.feature;

import junit.framework.AssertionFailedError;
import org.broad.igv.Globals;
import org.broad.igv.TestInformation;
import org.broad.igv.tools.IgvTools;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import javax.swing.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static junit.framework.Assert.*;

/**
 * User: jacob
 * Date: 2011/12/15
 */
public class FeatureDBTest {

    static String dataFileName = TestInformation.DATA_DIR + "/genomes/hg18.genome";
    public static final int LARGE = 500;

    private static final String CHECK_STR = "ABC";
    private static final int EXPECTED = 50;


    //Not a unit test
    public static void junk(String[] args) {
        //String CHECK_STR = CHECK_STR;
        //Map<String, NamedFeature> fMap = FeatureDB.getFeatures(CHECK_STR);

        String[] Options = {"Option1", "Option2", "Option3", "48549"};
        // Create the JList containing the items:
        JList ls = new JList(Options);
        ls.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

        // On the square in a scrollPane size desired:
        JScrollPane scls = new JScrollPane(ls);

        // It sets the initial val:
        ls.setSelectedValue(Options[0], true);

        // Create the content of our dialogue:
        // Message 1 + scrollpane containing the ls:
        Object[] message = {"Choose your contact", ls};

        // Use showOptionDialog () interface which offers the most complete
        int resp = JOptionPane.showConfirmDialog(
                null, message, "Contact List.",
                JOptionPane.OK_CANCEL_OPTION, JOptionPane.PLAIN_MESSAGE, null);

        // And we treat the return val:
        String val = null;
        if (resp == JOptionPane.OK_OPTION) {
            val = ls.getSelectedValue().toString();
        }
        //return val;

    }


    @Before
    public void setUp() throws Exception {
        Globals.setHeadless(true);
        IgvTools.loadGenome(dataFileName, true);
    }

    @After
    public void tearDown() {
        FeatureDB.clearFeatures();
    }

    @Test
    public void testFeaturesMap() throws Exception {
        Map<String, NamedFeature> fMap = FeatureDB.getFeaturesMap(CHECK_STR);

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

        Thread read = new Thread(new Runnable() {
            public void run() {
                try {
                    List<NamedFeature> features = FeatureDB.getFeaturesList(CHECK_STR, LARGE);
                    for (NamedFeature f : features) {
                        //Check for data corruption
                        assertTrue(f.getName().startsWith(CHECK_STR));
                    }
                    assertEquals(EXPECTED, features.size());
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

        List<NamedFeature> features = FeatureDB.getFeaturesList(CHECK_STR, LARGE);
        assertEquals(0, features.size());

        if (map.containsKey(0)) {
            AssertionFailedError e = map.get(0);
            throw e;
        }

    }


}
