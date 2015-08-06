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

import junit.framework.Assert;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.GenomeListItem;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.Track;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.fest.swing.fixture.FrameFixture;
import org.fest.swing.fixture.JButtonFixture;
import org.fest.swing.fixture.JComboBoxFixture;
import org.fest.swing.fixture.JPanelFixture;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Collection;
import java.util.List;

import static org.junit.Assert.*;

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

    /**
     * Test loading a UCSC session, with some files that don't exist and some that do
     *
     * @throws Exception
     */
    @Test
    public void testLoadUCSCSessionBadFiles() throws Exception {
        String sessionPath = TestUtils.DATA_DIR + "sessions/nonexistent_files.session";
        IGV igv = IGV.getInstance();

        TestUtils.loadSession(igv, sessionPath);

        Assert.assertEquals(4, igv.getVisibleTrackCount());
        int testBedTracks = 0;
        for (Track track : igv.getAllTracks()) {
            testBedTracks += track.getName().contains("Test Bed") ? 1 : 0;
        }
        Assert.assertEquals(2, testBedTracks);
    }

    @Ignore //Never seems to work on our testrunner for reasons unrelated to the code
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
        //Make sure frame has focus, or else homeButton won't work
        JButtonFixture homeButton = frameFixture.button("homeButton");

        homeButton.focus();
        homeButton.click();

        homeButton.focus();
        homeButton.click();

        igv.waitForNotify(500);

        Assert.assertEquals(Globals.CHR_ALL, frame.getChrName());

        //In all genome view these should be the same
        assertEquals(frame.getChromosomeLength(), frame.getCurrentRange().getEnd());
        Assert.assertEquals(0.0, frame.getOrigin());
    }

    @Test
    public void testLoadNewGenomeByPath() throws Exception {
        String sessionPath = TestUtils.DATA_DIR + "sessions/canFam2_local.xml";
        String genomeId = "canFam2.unittest";

        TestUtils.loadSession(IGV.getInstance(), sessionPath);
        assertEquals(genomeId, GenomeManager.getInstance().getGenomeId());
    }

    /**
     * Test loading a genome the user hasn't loaded before,
     * by id (available from server)
     *
     * @throws Exception
     */
    @Test
    public void testLoadNewGenomeById() throws Exception {
        Collection<GenomeListItem> currentGenomes = GenomeManager.getInstance().getGenomes();
        String genomeId = "canFam2";
        for (GenomeListItem genomeListItem : currentGenomes) {
            assertNotSame(genomeId, genomeListItem.getId());
        }
        String sessionPath = TestUtils.DATA_DIR + "sessions/canFam2_server.xml";
        TestUtils.loadSession(IGV.getInstance(), sessionPath);

        assertEquals(genomeId, GenomeManager.getInstance().getGenomeId());

        //CpG islands, RefSeq genes, ReferenceSequence (has height 0)
        assertEquals(3, IGV.getInstance().getVisibleTrackCount());
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
