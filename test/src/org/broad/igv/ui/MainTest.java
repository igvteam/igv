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

package org.broad.igv.ui;

import org.broad.igv.AbstractHeadedTest;
import org.broad.igv.PreferenceManager;
import org.broad.igv.track.Track;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.util.TestUtils;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

import static org.junit.Assert.*;

/**
 * Test of class main. In general will use this to see that IGV
 * starts properly, given startup parameters. Since we will be starting
 * IGV in a non-standard way, do NOT inherit from AbstractHeadedTest
 * User: jacob
 * Date: 2012/05/04
 */
public class MainTest {

    @Rule
    public TestRule testTimeout = new Timeout((int) 1e5);

    @BeforeClass
    public static void setUpClass() throws Exception {
        AbstractHeadedTest.assumeNotHeadless();
        TestUtils.setUpTestEnvironment();
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
        AbstractHeadedTest.tearDownClass();
    }

    /**
     * Test that loading IGV with a startup file and
     * locus loads that file and locus
     *
     * @throws Exception
     */
    @Test
    public void testFileLocusArgs() throws Exception {

        String trackName = "NA12878.SLX.bam";
        String filePath = "http://www.broadinstitute.org/igvdata/1KG/pilot2Bams/" + trackName;
        String genome = "hg18";
        String locus = "chr1:64,098,103-64,098,175";

        String[] args = new String[]{filePath, locus, "-g", genome};

        //Need to wait for IGV to load file, genome, and move to locus
        IGV igv = startMain(args, 60000);

        assertEquals(genome, igv.getGenomeManager().getGenomeId());
        boolean trackFound = false;
        for (Track track : igv.getAllTracks()) {
            trackFound |= track.getName().equals(trackName);
        }

        assertTrue("Expected " + trackName + " to be loaded, but it wasn't", trackFound);

        String actLocus = FrameManager.getDefaultFrame().getFormattedLocusString();
        assertEquals(actLocus, locus);
    }

    /**
     * Test loading a genome not in the display list
     *
     * @throws Exception
     */
    @Test
    public void testLoadGenome() throws Exception {

        String genomeId = "mm7";
        String[] genomes = PreferenceManager.getInstance().getGenomeIdDisplayList();
        for (String gen : genomes) {
            assertNotSame("Bad test setup, test genome in display list", gen, genomeId);
        }

        String[] args = new String[]{"-g", genomeId};
        IGV igv = startMain(args, 10000);

        assertEquals(igv.getGenomeManager().getGenomeId(), genomeId);
    }

    private IGV startMain(String[] args, long timeout) {
        Main.main(args);
        IGV igv = IGV.getInstance();
        assertTrue(igv.waitForNotify(timeout));
        return igv;
    }


}
