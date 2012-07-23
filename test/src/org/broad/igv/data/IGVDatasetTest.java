/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.data;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @author jrobinso
 */
public class IGVDatasetTest extends AbstractHeadlessTest {

    String cnFile = TestUtils.DATA_DIR + "igv/MIP_44.cn";

    /**
     * Test loading a dataset and accessing some random data cells
     */
    @Test
    public void testDataset() {


        IGVDataset ds = new IGVDataset(new ResourceLocator(cnFile), genome);

        // Get the start locations and data from sample yw280-4_44 on chr 1 
        int[] startLocations = ds.getStartLocations("chr1");
        float[] data = ds.getData("yw280-4_44", "chr1");
        assertEquals(startLocations.length, data.length);

        // Check the first data point
        assertEquals(794055, startLocations[0]);
        assertEquals(1.86971, data[0], 1.0e-4);

        // Check the last data point
        int end = startLocations.length - 1;
        assertEquals(245171540, startLocations[end]);
        assertEquals(3.7871208, data[end], 1.0e-4);

    }

    /**
     * Test of the parser getHeadings method
     */
    @Test
    public void testGetHeadings() {
        String firstHeading = "yw280-4_44";
        String secondHeading = "yw280-4_48";
        String lastHeading = "DN30_3_03";
        String headingsLine = "SNP	Chromosome	PhysicalPosition	yw280-4_44	yw280-4_48	yw280-5_40	yw280-5_41	yw280-5_42	yw280-5_43	yw280-5_44	yw280-5_46	yw280-5_47	yw280-5_48	yw280-5_45	DN30_01	DN30_02	DN30_03	DN30_04	DN30_06	DN30_07	DN30_09	DN30_10	DN30_11	DN30_13	DN30_14	DN30_15	DN30_16	DN30_17	DN30_18	DN30_19	DN30_20	DN30_21	DN30_32	DN30_33	DN30_34	DN30_35	DN30_36	DN30_41	DN30_42	DN30_43	DN30_44	DN30_45	DN30_46	DN30_47	DN30_48	DN30_3_02	DN30_3_03";

        String[] tokens = headingsLine.split("\t");

        IGVDataset ds = new IGVDataset(new ResourceLocator(cnFile), genome);
        IGVDatasetParser parser = new IGVDatasetParser(new ResourceLocator(cnFile), genome);

        //parser.scan necessary to set values needed in getHeadings
        List<ChromosomeSummary> summaries = parser.scan(ds);
        String[] headings = parser.getHeadings(tokens, 1);

        assertEquals(firstHeading, headings[0]);
        assertEquals(secondHeading, headings[1]);
        assertEquals(lastHeading, headings[headings.length - 1]);
    }

    /**
     * Test of scan method. Test that we find
     * all chromosomes we expect
     */
    @Test
    public void testScanDataset() {
        IGVDataset ds = new IGVDataset(new ResourceLocator(cnFile), genome);
        IGVDatasetParser parser = new IGVDatasetParser(new ResourceLocator(cnFile), genome);
        List<ChromosomeSummary> summaries = parser.scan(ds);

        assertEquals(24, summaries.size());
        Set<String> chromos = new TreeSet<String>();
        ArrayList<String> chromos_L = new ArrayList<String>(24);
        for (int ii = 1; ii <= 22; ii++) {
            chromos_L.add("chr" + ii);
        }
        chromos_L.add("chrX");
        chromos_L.add("chrY");
        chromos.addAll(chromos_L);

        for (ChromosomeSummary cSum : summaries) {
            assertTrue(chromos.contains(cSum.getName()));
        }

    }


}
