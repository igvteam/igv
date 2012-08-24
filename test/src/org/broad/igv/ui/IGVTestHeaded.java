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

import junit.framework.Assert;
import org.broad.igv.AbstractHeadedTest;
import org.broad.igv.Globals;
import org.broad.igv.track.Track;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.fest.swing.fixture.FrameFixture;
import org.fest.swing.fixture.JComboBoxFixture;
import org.fest.swing.fixture.JPanelFixture;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012/02/08
 */
public class IGVTestHeaded extends AbstractHeadedTest {

    @Test
    public void testLoadSessionBatch() throws Exception {
        try {
            Globals.setBatch(true);
            tstLoadSession();
        } finally {
            Globals.setBatch(false);
        }

    }

    @Test
    public void testLoadSessionNoBatch() throws Exception {
        Globals.setBatch(false);
        tstLoadSession();
    }

    public void tstLoadSession() throws Exception {
        //Pretty basic, but at some point loading this view
        //gave a class cast exception
        String sessionPath = TestUtils.DATA_DIR + "sessions/CCLE_testSession_chr2.xml";
        IGV igv = IGV.getInstance();

        TestUtils.loadSession(igv, sessionPath);

        Assert.assertEquals("chr2", FrameManager.getDefaultFrame().getChrName());
        Assert.assertEquals(1, FrameManager.getDefaultFrame().getCurrentRange().getStart());

        int rangeDiff = Math.abs(FrameManager.getDefaultFrame().getChromosomeLength() - FrameManager.getDefaultFrame().getCurrentRange().getEnd());
        assertTrue(rangeDiff < 3);

        Assert.assertEquals(1461, igv.getAllTracks().size());
    }

    @Test
    public void testHome() throws Exception {
        IGV igv = IGV.getInstance();
        ReferenceFrame frame = FrameManager.getDefaultFrame();
        String chr = "chr1";
        int start = 5;
        int end = 5000;
        int limit = 2;

        frame.jumpTo(chr, start, end);

        Assert.assertEquals(chr, frame.getChrName());
        assertTrue(Math.abs(frame.getCurrentRange().getStart() - start) < limit);
        assertTrue(Math.abs(frame.getCurrentRange().getEnd() - end) < limit);


        FrameFixture frameFixture = new FrameFixture(IGV.getMainFrame());
        frameFixture.button("homeButton").click();

        IGV.getInstance().waitForNotify(500);

        Assert.assertEquals(Globals.CHR_ALL, frame.getChrName());

        //In all genome view these should be the same
        assertEquals(frame.getChromosomeLength(), frame.getCurrentRange().getEnd());
        Assert.assertEquals(0.0, frame.getOrigin());


    }

    /**
     * Basic test showing usage of FEST and checking combo box
     *
     * @throws Exception
     */
    //@Test
    public void scratchTestFEST() throws Exception {

        FrameFixture frame = new FrameFixture(IGV.getMainFrame());
        JPanelFixture contentFixture = frame.panel("contentPane");

        JPanelFixture commandBar = frame.panel("igvCommandBar");
        JComboBoxFixture chromoBox = frame.comboBox("chromosomeComboBox");

        String[] chromos = commandBar.comboBox("chromosomeComboBox").contents();
        Assert.assertEquals(26, chromos.length);
    }


    public void testExomeView() throws Exception {
        String file = "http://www.broadinstitute.org/igvdata/1KG/pilot2Bams/NA12891.SLX.bam";
        List<Track> tracks = IGV.getInstance().load(new ResourceLocator(file));
        Thread.sleep(10000);
        //TestUtils.loadSession(igv, TestUtils.DATA_DIR + "sessions/slx_ceu_father.xml");

        Assert.assertEquals(2, tracks.size());

        FrameManager.setExomeMode(true, true);
        IGV.getInstance().resetFrames();

        String locus = "chr7:55,208,260-55,240,460";
        igv.goToLocus(locus);

        //File out = new File(TestUtils.DATA_DIR, "testsnap.png");
        //SnapshotUtilities.doComponentSnapshot(IGV.getMainFrame(), out, SnapshotFileChooser.SnapshotFileType.PNG);


    }
}
