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

package org.broad.igv.feature;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeDescriptor;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.IGV;
import org.junit.AfterClass;
import static org.junit.Assert.assertEquals;
import org.junit.BeforeClass;
import org.junit.Test;

import java.util.Arrays;

/**
 * @author jrobinso
 */
public class GenomeManagerTest {

    static GenomeManager genomeManager;

    public GenomeManagerTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        genomeManager = IGV.getInstance().getGenomeManager();
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }



    @Test
    public void testSortChromosomes() {

        String [] chrs = {"chr12", "chr10", "chrMT",  "chr1", "chrLongName", "chrLongName1"};
        String [] expectedResult = {"chr1", "chr10", "chr12", "chrLongName1", "chrLongName", "chrMT"};

        Arrays.sort(chrs, new Genome.ChromosomeComparator());
        for(int i=0; i<chrs.length; i++) {
            assertEquals(expectedResult[i], chrs[i]);
        }
        System.out.println();

        chrs = new String [] {"scaffold_v2_10414", "scaffold_v2_100", "scaffold_v2_101", "scaffold_v2_10415"};
        expectedResult = new String [] {"scaffold_v2_100", "scaffold_v2_101", "scaffold_v2_10414", "scaffold_v2_10415"};
        Arrays.sort(chrs, new Genome.ChromosomeComparator());
        for(int i=0; i<chrs.length; i++) {
            assertEquals(expectedResult[i], chrs[i]);
        }
    }

    @Test
    public void testGetDefaultGenomeDescriptor() {
        GenomeDescriptor gd =  genomeManager.getDefaultGenomeDescriptor();
        assertEquals("hg18", gd.getId());
       
    }


}