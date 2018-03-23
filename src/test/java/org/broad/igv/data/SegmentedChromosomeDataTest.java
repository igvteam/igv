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
