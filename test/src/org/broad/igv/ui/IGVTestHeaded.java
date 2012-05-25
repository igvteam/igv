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
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.RegionScoreType;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.igv.util.Utilities;
import org.junit.*;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012/02/08
 */
public class IGVTestHeaded extends AbstractHeadedTest{

    @Test
    public void testLoadSession() throws Exception{
        //Pretty basic, but at some point loading this view
        //gave a class cast exception
        String sessionPath = TestUtils.DATA_DIR + "sessions/CCLE_testSession_chr2.xml";
        IGV igv = IGV.getInstance();

        TestUtils.loadSession(igv, sessionPath);

        assertEquals("chr2", FrameManager.getDefaultFrame().getChrName());
        assertEquals(1, FrameManager.getDefaultFrame().getCurrentRange().getStart());

        int rangeDiff = Math.abs(FrameManager.getDefaultFrame().getChromosomeLength() - FrameManager.getDefaultFrame().getCurrentRange().getEnd());
        assertTrue(rangeDiff < 3);
    }

}
