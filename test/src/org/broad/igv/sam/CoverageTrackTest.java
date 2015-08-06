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
        int startingNumTracks = igv.getAllTracks().size();
        String bamPath = TestUtils.LARGE_DATA_DIR + "HG00171.hg18.bam";
        Collection<ResourceLocator> locators = Arrays.asList(new ResourceLocator(bamPath));
        igv.loadResources(locators);
        List<Track> tracks = igv.getAllTracks();
        assertEquals(startingNumTracks + 2, tracks.size());

        String startLocString = "chr1:151666494-151666594";

        String midLocString = "chr1:153148479-153148579";
        Locus midLocus = Locus.fromString(midLocString);

        String destLocString = "chr1:155232055-155232155";
        Locus destLocus = Locus.fromString(destLocString);


        AlignmentTrack alTrack = (AlignmentTrack) tracks.get(1);
        AlignmentDataManager dataManager = (alTrack).getDataManager();

        igv.goToLocus(startLocString);

        checkLocus(dataManager, destLocus, false);
        checkLocus(dataManager, midLocus, false);

        igv.removeTracks(Arrays.<Track>asList(alTrack));
        igv.goToLocus(destLocString);
        Thread.sleep(3000);

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
