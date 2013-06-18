/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.ui;

import org.broad.igv.PreferenceManager;
import org.broad.igv.sam.AlignmentBlock;
import org.broad.igv.sam.AlignmentDataManager;
import org.broad.igv.sam.AlignmentInterval;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.track.Track;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.*;

/**
 * @author jacob
 * @date 2013-Jun-18
 */
public class PreferencesEditorTest extends AbstractHeadedTest {

    /**
     * Certain preferences should force a reload when they change,
     * mainly those associated with Alignments
     *
     * @throws Exception
     */
    @Test
    public void testChangeAlignmentReload() throws Exception {

        String showSC = PreferenceManager.SAM_SHOW_SOFT_CLIPPED;

        PreferenceManager.getInstance().put(showSC, "false");
        assertFalse(PreferenceManager.getInstance().getAsBoolean(showSC));

        //Load a data file which has a soft clipped block
        String path = TestUtils.DATA_DIR + "sam/hardSoftClip.sam";
        TestUtils.createIndex(path);
        List<Track> tracks = IGV.getInstance().load(new ResourceLocator(path));
        AlignmentTrack alTrack = (AlignmentTrack) tracks.get(1);

        //NPE happens when event is broadcast, clutters up stacktrace
        alTrack.getDataManager().getEventBus().unregister(tracks.get(0));

        AlignmentBlock[] blocks0 = getBlocks(alTrack, "chr1", 59300, 59400);
        assertEquals(1, blocks0.length);
        assertFalse(blocks0[0].isSoftClipped());

        PreferencesEditor dialog = new PreferencesEditor(IGV.getMainFrame(), true);

        dialog.updatedPreferenceMap.put(showSC, "true");
        PreferenceManager.getInstance().put(showSC, "true");
        dialog.okButton.doClick();

        assertTrue(PreferenceManager.getInstance().getAsBoolean(showSC));

        AlignmentBlock[] blocks1 = getBlocks(alTrack, "chr1", 59300, 59400);
        assertEquals(2, blocks1.length);
        assertTrue(blocks1[0].isSoftClipped());
        assertFalse(blocks1[1].isSoftClipped());

    }

    /**
     * Loads the AlignmentBlocks from the specified range. Cached value is used if available,
     * otherwise data is loaded
     *
     * @param alTrack
     * @param chr
     * @param start
     * @param end
     * @return
     */
    private AlignmentBlock[] getBlocks(AlignmentTrack alTrack, String chr, int start, int end) {

        AlignmentDataManager dataManager = alTrack.getDataManager();
        AlignmentInterval interval = dataManager.getLoadedInterval(FrameManager.getDefaultFrame().getName());
        if (interval == null) {
            dataManager.loadAlignments(chr, start, end, null, null);
            interval = dataManager.getLoadedInterval(FrameManager.getDefaultFrame().getName());
        }

        return interval.getAlignmentIterator().next().getAlignmentBlocks();

    }
}
