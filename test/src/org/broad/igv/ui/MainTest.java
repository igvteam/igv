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

import org.broad.igv.track.Track;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.util.TestUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012/05/04
 */
public class MainTest {

    @Before
    public void setUp() throws Exception {

    }

    @After
    public void tearDown() throws Exception {

    }

    /**
     * Test that loading IGV with a startup file and
     * locus loads that file and locus
     *
     * @throws Exception
     */
    @Test
    public void testFileLocusArgs() throws Exception {
        TestUtils.assumeNotHeadless();

        String trackName = "NA12878.SLX.bam";
        String filePath = "http://www.broadinstitute.org/igvdata/1KG/pilot2Bams/" + trackName;
        String genome = "hg18";
        String locus = "chr1:64,098,103-64,098,175";

        String[] args = new String[]{filePath, locus, "-g", genome};

        Main.main(args);

        IGV igv = IGV.getInstance();

        //Need to wait for IGV to load file, genome,  and move to locus
        assertTrue(igv.waitForNotify(60000));

        assertEquals(genome, igv.getGenomeManager().getGenomeId());
        boolean trackFound = false;
        for (Track track : igv.getAllTracks(false)) {
            trackFound |= track.getName().equals(trackName);
        }

        assertTrue("Expected " + trackName + " to be loaded, but it wasn't", trackFound);

        String actLocus = FrameManager.getDefaultFrame().getFormattedLocusString();
        assertEquals(actLocus, locus);


    }


}
