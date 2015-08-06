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
        AlignmentInterval interval = dataManager.getLoadedInterval(FrameManager.getDefaultFrame().getCurrentRange());
        if (interval == null) {
            dataManager.loadAlignments(chr, start, end, null, null);
            interval = dataManager.getLoadedInterval(FrameManager.getDefaultFrame().getCurrentRange());
        }

        return interval.getAlignmentIterator().next().getAlignmentBlocks();

    }
}
