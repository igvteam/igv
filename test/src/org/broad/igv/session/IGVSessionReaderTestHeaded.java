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

package org.broad.igv.session;

import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.sam.CoverageTrack;
import org.broad.igv.track.Track;
import org.broad.igv.ui.AbstractHeadedTest;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.lang.reflect.Method;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/** Test reading sessions, make sure attributes get parsed properly
 *
 * TODO USE LOCAL DATA FILES
 * User: jacob
 * Date: 2013-Jan-07
 */
public class IGVSessionReaderTestHeaded extends AbstractHeadedTest{

    @Test
    public void testLoadCoverageTrackSession() throws Exception{
        String path = TestUtils.DATA_DIR + "sessions/coverage_snpThreshold.xml";
        IGV.getInstance().doRestoreSession(path, null, false);

        //We should have 1 coverage track
        CoverageTrack track = null;
        for(Track tk: IGV.getInstance().getAllTracks()){
            if(tk instanceof CoverageTrack){
                track = (CoverageTrack) tk;
                break;
            }

        }

        assertTrue(track.isShowReference());
        assertEquals(0.5, track.getSnpThreshold(), 1e-5);
        assertTrue(track.isAutoScale());
    }

    /**
     * Loads path, returns first alignment track
     * @param path
     * @return
     */
    private AlignmentTrack getAlignmentTrack(String path){
        IGV.getInstance().doRestoreSession(path, null, false);

        AlignmentTrack alTrack = null;
        for(Track tk: IGV.getInstance().getAllTracks()){
            if(tk instanceof AlignmentTrack){
                alTrack = (AlignmentTrack) tk;
            }
        }

        return alTrack;
    }

    @Test
    public void testLoadAlignmentTrackSession() throws Exception{
        String path = TestUtils.DATA_DIR + "sessions/HG00171_wlocus_v5.xml";
        AlignmentTrack alTrack = getAlignmentTrack(path);

        CoverageTrack covTrack = alTrack.getCoverageTrack();
        assertTrue(alTrack.isShowSpliceJunctions());

        assertEquals(CoverageTrack.DEFAULT_SHOW_REFERENCE, covTrack.isShowReference());
        assertEquals(CoverageTrack.DEFAULT_AUTOSCALE, covTrack.isAutoScale());
    }


    @Test
    public void testLoadRenderOptions_v4() throws Exception{
        String path = TestUtils.DATA_DIR + "sessions/HG00171_wlocus_v4.xml";
        tstLoadRenderOptions(path);
    }

    @Test
    public void testLoadRenderOptions_v5() throws Exception{
        String path = TestUtils.DATA_DIR + "sessions/HG00171_wlocus_v5.xml";
        tstLoadRenderOptions(path);
    }

    private void tstLoadRenderOptions(String path) throws Exception{
        AlignmentTrack alTrack = getAlignmentTrack(path);

        AlignmentTrack.RenderOptions ro = getRenderOptions(alTrack);

        assertEquals(AlignmentTrack.GroupOption.FIRST_OF_PAIR_STRAND, ro.getGroupByOption());
        assertEquals(AlignmentTrack.ColorOption.INSERT_SIZE, ro.getColorOption());
        assertEquals(55, ro.getMinInsertSize());
        assertEquals(1001, ro.getMaxInsertSize());
        assertEquals("", ro.getColorByTag());
    }

    /**
     * Get render options of the track, reflectively so we don't need to alter access controls
     * @param track
     * @return
     * @throws Exception
     */
    private AlignmentTrack.RenderOptions getRenderOptions(Track track) throws Exception{
        Method getRenderOptions = track.getClass().getDeclaredMethod("getRenderOptions");
        getRenderOptions.setAccessible(true);
        Object result = getRenderOptions.invoke(track);
        return (AlignmentTrack.RenderOptions) result;
    }
}
