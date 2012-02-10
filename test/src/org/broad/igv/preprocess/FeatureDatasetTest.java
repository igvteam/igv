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
package org.broad.igv.preprocess;

import org.broad.igv.data.expression.ExpressionDataset;
import org.broad.igv.data.expression.ExpressionFileParser;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.TestUtils;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;

import static org.junit.Assert.assertEquals;

/**
 * @author jrobinso
 */
public class FeatureDatasetTest {

    static String file = TestUtils.DATA_DIR + "/gct/affy_human_mod.gct";
    static ExpressionDataset dataset;

    public FeatureDatasetTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        TestUtils.setUpHeadless();
        Genome genome = TestUtils.loadGenome();
        ExpressionFileParser parser = new ExpressionFileParser(new File(file), null, genome);

        dataset = parser.createDataset();
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
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
