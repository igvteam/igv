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

package org.broad.igv.sam;

import org.broad.igv.Globals;
import org.broad.igv.feature.Locus;
import org.broad.igv.track.Track;
import org.broad.igv.ui.AbstractHeadedTest;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import static junit.framework.Assert.assertEquals;

/**
 * User: jacob
 * Date: 2012-Jul-23
 */
public class CoverageTrackTest extends AbstractHeadedTest {

    @Before
    public void setUp() throws Exception {
        super.setUp();
        Globals.setBatch(true);
    }

    @After
    public void tearDown() throws Exception {
        super.tearDown();
        Globals.setBatch(true);
    }


    /**
     * The purpose of this test is to see if the coverage track
     * still renders if a BAM file is loaded and then removed.
     * <p/>
     * We actually look at the dataManager, because there's
     * no good way to check what was rendered
     *
     * @throws Exception
     */
    @Test
    public void testAfterBAMRemoved() throws Exception {
        String bamPath = TestUtils.LARGE_DATA_DIR + "HG00171.hg18.bam";
        Collection<ResourceLocator> locators = Arrays.asList(new ResourceLocator(bamPath));
        igv.loadResources(locators);
        List<Track> tracks = igv.getAllTracks();
        assertEquals(2, tracks.size());

        String startLocString = "chr1:151666494-151666594";

        String midLocString = "chr1:153148479-153148579";
        Locus midLocus = new Locus(midLocString);

        String destLocString = "chr1:155232055-155232155";
        Locus destLocus = new Locus(destLocString);


        AlignmentTrack alTrack = (AlignmentTrack) tracks.get(1);
        AlignmentDataManager dataManager = (alTrack).getDataManager();

        igv.goToLocus(startLocString);

        checkLocus(dataManager, destLocus, false);
        checkLocus(dataManager, midLocus, false);

        igv.removeTracks(Arrays.<Track>asList(alTrack));
        igv.goToLocus(destLocString);

        checkLocus(dataManager, destLocus, true);

        //Check something that shouldn't be loaded, make
        //sure it doesn't get loaded when BOTH tracks removed
        igv.removeTracks(tracks);
        igv.goToLocus(midLocString);
        checkLocus(dataManager, midLocus, false);

    }

    private void checkLocus(AlignmentDataManager dataManager, Locus locus, boolean expectContain) {
        boolean contains = false;

        for (AlignmentInterval interval : dataManager.getLoadedIntervals()) {
            contains |= interval.contains(locus.getChr(), locus.getStart(), locus.getEnd());
        }

        assertEquals(expectContain, contains);
    }
}
