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
package org.broad.igv.preprocess;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.data.expression.ExpressionDataset;
import org.broad.igv.data.expression.ExpressionFileParser;
import org.broad.igv.util.TestUtils;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;

import static org.junit.Assert.assertEquals;

/**
 * @author jrobinso
 */
public class FeatureDatasetTest extends AbstractHeadlessTest {

    static String file = TestUtils.DATA_DIR + "gct/affy_human_mod.gct";
    static ExpressionDataset dataset;

    public FeatureDatasetTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        AbstractHeadlessTest.setUpClass();
        ExpressionFileParser parser = new ExpressionFileParser(new File(file), null, genome);

        dataset = parser.createDataset();
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
        AbstractHeadlessTest.tearDownClass();
        dataset = null;
    }

    /**
     * Test of getDataHeadings method, of class ExpressionDataset.
     * Excludes probe" and description columns
     */
    @Test
    public void getDataHeadings() {

        int expectedSize = 7;
        String[] headings = dataset.getTrackNames();
        assertEquals(expectedSize, headings.length);
    }


    /**
     * Test of getChromosomes method, of class ExpressionDataset.
     */
    @Test
    public void getChromosomes() {
        int nChromosomes = 24;
        String[] chromosomes = dataset.getChromosomes();
        assertEquals(nChromosomes, chromosomes.length);
    }

    /**
     * Test of getStartLocations method, of class ExpressionDataset.
     */
    @Test
    public void getStartLocations() {
        int expectedSize = 10;
        int[] locations = dataset.getStartLocations("chr1");
        assertEquals(expectedSize, locations.length);

        int firstLoc = 19883202;
        assertEquals(firstLoc, locations[0]);

        int lastLoc = 241359291;
        assertEquals(lastLoc, locations[expectedSize - 1]);

    }

    /**
     * Test of getEndLocations method, of class ExpressionDataset.
     */
    @Test
    public void getEndLocations() {
        int expectedSize = 10;
        int[] locations = dataset.getEndLocations("chr1");
        assertEquals(expectedSize, locations.length);

        int firstLoc = 19883459;
        assertEquals(firstLoc, locations[0]);

        int lastLoc = 241359329;
        assertEquals(lastLoc, locations[expectedSize - 1]);

    }

    /**
     * Test of getData method, of class ExpressionDataset.
     */
    @Test
    public void getData() {

        String chr = "chr10";
        String heading = "AML_1";

        int expectedSize = 3;
        float[] data = dataset.getData(heading, chr);
        assertEquals(expectedSize, data.length);

        float value = 573.0001f;
        assertEquals(value, data[0], 0.0000001f);

        float lastValue = 573.0001f;
        assertEquals(lastValue, data[expectedSize - 1], 0.0000001f);


    }
}
