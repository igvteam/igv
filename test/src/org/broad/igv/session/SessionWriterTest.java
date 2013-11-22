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

package org.broad.igv.session;

import org.broad.igv.track.DataTrack;
import org.broad.igv.track.MergedTracks;
import org.broad.igv.ui.AbstractHeadedTest;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.action.SaveSessionMenuAction;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.File;
import java.util.List;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * @author jacob
 * @date 2013-Nov-22
 */
public class SessionWriterTest extends AbstractHeadedTest {

    @Test
    public void testSaveOverlaySession() throws Exception{
        String sessionpath = TestUtils.DATA_DIR + "sessions/hind_gistic_overlay.xml";
        rewriteRestoreSession(sessionpath);

        File targetFile = new File(TestUtils.TMP_OUTPUT_DIR, "overlay_session.xml");
        SaveSessionMenuAction.saveSession(IGV.getInstance(), targetFile);

        rewriteRestoreSession(targetFile.getAbsolutePath());


        List<DataTrack> dataTracks = IGV.getInstance().getDataTracks();
        assertEquals(3, dataTracks.size());
        assertTrue(dataTracks.get(0) instanceof MergedTracks);
    }
}
