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
package org.broad.igv.data;

import org.junit.*;
import static org.junit.Assert.assertEquals;
import org.broad.igv.data.seg.SegmentedChromosomeData;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

/**
 * @author jrobinso
 */
public class SegmentedChromosomeDataTest {

    SegmentedChromosomeData testData;

    public SegmentedChromosomeDataTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
        String[] samples = {"sample 1", "sample 2", "sample 3"};

        Map<String, int[]> starts = new HashMap();
        Map<String, int[]> ends = new HashMap();
        Map<String, float[]> values = new HashMap();
        for (String s : samples) {
            starts.put(s, new int[]{1, 50, 100});
            ends.put(s, new int[]{25, 75, 125});
            values.put(s, new float[]{1.0f, 2.0f, 3.0f});
        }

        testData = new SegmentedChromosomeData(samples, starts, ends, values);
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of serialize method, of class SegmentedChromosomeData.
     */
    @Test
    public void serialize() throws Exception {

        File testFile = File.createTempFile("test", "bin");
        testFile.deleteOnExit();


        OutputStream stream = new FileOutputStream(testFile);
        testData.serialize(stream);
        stream.flush();
        stream.close();

        SegmentedChromosomeData copy = new SegmentedChromosomeData();
        InputStream is = new FileInputStream(testFile);
        copy.deserialize(is);
        is.close();

        String[] expectedSamples = testData.getSampleNames();
        String[] samples = copy.getSampleNames();
        assertEquals(expectedSamples.length, samples.length);

        for (int i = 0; i < samples.length; i++) {
            assertEquals(expectedSamples[i], samples[i]);
        }

        for (String s : samples) {
            int[] expectedStarts = testData.getStartLocations(s);
            int[] starts = copy.getStartLocations(s);
            assertEquals(expectedStarts.length, starts.length);
            for (int i = 0; i < starts.length; i++) {
                assertEquals(expectedStarts[i], starts[i]);
            }

            int[] expectedEnds = testData.getEndLocations(s);
            int[] ends = copy.getEndLocations(s);
            assertEquals(expectedEnds.length, ends.length);
            for (int i = 0; i < starts.length; i++) {
                assertEquals(expectedEnds[i], ends[i]);
            }

            float[] expectedValues = testData.getValues(s);
            float[] values = copy.getValues(s);
            assertEquals(values.length, values.length);
            for (int i = 0; i < starts.length; i++) {
                assertEquals(expectedValues[i], values[i], 1.0e-6);
            }
        }


    }
}
