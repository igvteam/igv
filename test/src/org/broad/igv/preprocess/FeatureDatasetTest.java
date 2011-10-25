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
import org.broad.igv.feature.genome.GenomeManager;
import org.junit.AfterClass;
import static org.junit.Assert.assertEquals;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;

/**
 * @author jrobinso
 */
public class FeatureDatasetTest {

    static String file = "/Volumes/xchip_tcga/gbm/visualization/testfiles/expression/WETime.rma.mapped";
    static ExpressionDataset dataset;
    static String genomeId = "mm8";

    public FeatureDatasetTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        Genome genome =   (new GenomeManager()).loadGenome("", null);
        ExpressionFileParser parser = new ExpressionFileParser(new File(file), null, genome);

        dataset = parser.createDataset();


    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    /**
     * Test of getDataHeadings method, of class ExpressionDataset.
     */
    @Test
    public void getDataHeadings() {
        int expectedSize = 16;
        String[] headings = dataset.getTrackNames();
        assertEquals(expectedSize, headings.length);
    }


    /**
     * Test of getChromosomes method, of class ExpressionDataset.
     */
    @Test
    public void getChromosomes() {
        int nChromosomes = 22;
        String[] chromosomes = dataset.getChromosomes();
        assertEquals(nChromosomes, chromosomes.length);
    }

    /**
     * Test of getStartLocations method, of class ExpressionDataset.
     */
    @Test
    public void getStartLocations() {
        int expectedSize = 2698;
        int[] locations = dataset.getStartLocations("chr1");
        assertEquals(expectedSize, locations.length);

        int firstLoc = 3063333;
        assertEquals(firstLoc, locations[0]);

        int lastLoc = 196841902;
        assertEquals(lastLoc, locations[expectedSize - 1]);

    }

    /**
     * Test of getEndLocations method, of class ExpressionDataset.
     */
    @Test
    public void getEndLocations() {
        int expectedSize = 2698;
        int[] locations = dataset.getEndLocations("chr1");
        assertEquals(expectedSize, locations.length);

        int firstLoc = 3064403;
        assertEquals(firstLoc, locations[0]);

        int lastLoc = 196877438;
        assertEquals(lastLoc, locations[expectedSize - 1]);

    }

    /**
     * Test of getData method, of class ExpressionDataset.
     */
    @Test
    public void getData() {

        String chr = "chr1";
        String heading = "SG20080102.B02.CEL";

        int expectedSize = 2698;
        float[] data = dataset.getData(heading, chr);
        assertEquals(expectedSize, data.length);

        float value = 3.964889866f;
        assertEquals(value, data[0], 0.0000001f);

        float lastValue = 3.219433646f;
        assertEquals(lastValue, data[expectedSize - 1], 0.0000001f);


    }
}
