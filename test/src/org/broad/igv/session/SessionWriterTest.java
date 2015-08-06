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
