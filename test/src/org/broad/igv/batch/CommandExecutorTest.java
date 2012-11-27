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

package org.broad.igv.batch;

import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.RegionScoreType;
import org.broad.igv.track.Track;
import org.broad.igv.ui.AbstractHeadedTest;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.IGVTestHeadless;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.util.TestUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Timer;

import static junit.framework.Assert.*;

/**
 * User: jacob
 * Date: 2012/03/21
 */
public class CommandExecutorTest extends AbstractHeadedTest {

    CommandExecutor exec = new CommandExecutor();
    private final String snapshotDir = TestUtils.TMP_OUTPUT_DIR;

    @Rule
    public TestRule testTimeout = new Timeout((int) 180000);

    @Before
    public void setUp() throws Exception {
        super.setUp();
        Globals.setBatch(true);
        igv.loadGenome(TestUtils.defaultGenome, null);
        igv.removeTracks(igv.getAllTracks());
        igv.getSession().clearRegionsOfInterest();
        exec.setSnapshotDirectory(snapshotDir);
    }

    @After
    public void tearDown() throws Exception {
        igv.removeTracks(igv.getAllTracks());
        Globals.setBatch(false);
    }

    @Test
    public void testRegion() throws Exception {

        String regionStr = "chr1:50-1000";
        exec.execute("region " + regionStr);

        Collection<RegionOfInterest> regions = IGV.getInstance().getSession().getAllRegionsOfInterest();
        assertEquals(1, regions.size());

        ArrayList<RegionOfInterest> regionsAL = (ArrayList<RegionOfInterest>) regions;
        RegionOfInterest region = regionsAL.get(0);
        assertEquals("chr1", region.getChr());
        assertEquals(49, region.getStart());
        assertEquals(1000, region.getEnd());

    }

    @Test
    public void testSetMaxDepth() throws Exception {
        setCheckMaxDepth(5);
        setCheckMaxDepth(50);
    }

    private void setCheckMaxDepth(int maxDepth) {
        String res = exec.execute("maxDepth " + maxDepth);
        assertFalse(res.contains("ERROR"));
        int newMaxDepth = PreferenceManager.getInstance().getAsInt(PreferenceManager.SAM_SAMPLING_COUNT);
        assertEquals(maxDepth, newMaxDepth);

    }

    @Test
    public void testSortByRegionScoreType() throws Exception {
        Timer deadlockChecker = TestUtils.startDeadlockChecker(1000);
        String sessionPath = TestUtils.DATA_DIR + "sessions/BRCA_loh2.xml";
        TestUtils.loadSession(igv, sessionPath);
        Collection<RegionOfInterest> rois = igv.getSession().getAllRegionsOfInterest();

        List<Track> tracks;
        int count = 0;
        for (RegionOfInterest roi : rois) {
            for (RegionScoreType type : RegionScoreType.values()) {
                igv.sortAllTracksByAttributes(new String[]{"NAME"}, new boolean[]{false});
                String typeStr = type.toString().toUpperCase();
                if (count % 2 == 0) {
                    typeStr = typeStr.toLowerCase();
                }
                String resp = exec.execute("sort " + typeStr + " " + roi.getLocusString());
                assertEquals("OK", resp);

                tracks = igv.getAllTracks();
                IGVTestHeadless.checkIsSorted(tracks, roi, type, FrameManager.getDefaultFrame().getZoom());
                count++;
            }
        }
        deadlockChecker.cancel();
        deadlockChecker.purge();
    }

    private final String outFileBase = "testSnap";

    @Test
    public void testSnapShotPng() throws Exception {
        String outFileName = outFileBase + ".Png";
        tstSnapshot(outFileName);
    }

    @Test
    public void testSnapshotTracksOnly() throws Exception {
        String outFileName = outFileBase + "_track.png";
        tstSnapshot(outFileName, true, "trackpanels");
    }

    @Test
    public void testSnapShotJpeg() throws Exception {
        tstSnapshot(outFileBase + ".jpeg");
    }

    @Test
    public void testSnapShotJpg() throws Exception {
        tstSnapshot(outFileBase + ".jpg");
    }

    @Test
    public void testSnapShotSvg() throws Exception {
        String outFileName = outFileBase + ".svG";
        tstSnapshot(outFileName);
    }

    @Test
    public void testSnapShotFails() throws Exception {
        String[] exts = new String[]{"abc", "svt", "pnq"};
        for (String ext : exts) {
            String outFileName = outFileBase + "." + ext;
            tstSnapshot(outFileName, false, null);
        }
    }


    public File tstSnapshot(String outFileName) throws Exception {
        return tstSnapshot(outFileName, true, null);
    }

    public File tstSnapshot(String outFileName, boolean shouldSucceed, String moreargs) throws Exception {
        File out = new File(snapshotDir, outFileName);
        assertFalse(out.exists());

        String toexec = "snapshot " + outFileName;
        if (moreargs != null && moreargs.length() > 0) {
            toexec += " " + moreargs;
        }
        exec.execute(toexec);

        assertEquals(shouldSucceed, out.exists());

        return out;
    }

    @Test
    public void testLoadURL() throws Exception {
        //This is mostly to ruggedize test setup. The setup may load
        //reference/sequence tracks, we'd like to be able to change
        //test setup and not worry about this test.
        int beginTracks = igv.getAllTracks().size();

        String urlPath = "ftp://ftp.broadinstitute.org/distribution/igv/TEST/cpgIslands%20with%20spaces.hg18.bed";
        exec.loadFiles(urlPath, null, true, "hasSpaces");

        String localPath = TestUtils.DATA_DIR + "bed/test.bed";
        exec.loadFiles(localPath, null, true, null);

        assertEquals(2, igv.getAllTracks().size() - beginTracks);
    }

    @Test
    public void testSetDataRange() throws Exception {
        String dataFile = TestUtils.DATA_DIR + "igv/recombRate.ens.igv.txt";
        exec.loadFiles(dataFile, null, true, null);

        String[] goodArgSet = new String[]{"0,5.0 ", "0,1,5", "-1,0,1", "-1.32,10.21"};
        for (String arg : goodArgSet) {
            String resp = exec.execute("setDataRange " + arg);
            assertEquals("OK", resp);
        }

        String[] badArgSet = new String[]{"0 ", "-1,0,2,3", "o,1o"};
        for (String arg : badArgSet) {
            String resp = exec.execute("setDataRange " + arg);
            assertTrue(resp.toLowerCase().startsWith("error"));
        }


    }

    @Test
    public void testLoadGenomesById() throws Exception {
        String[] genomeIds = new String[]{"hg19", "mm10", "rn5", "canFam2", "bosTau7", "sacCer3", "WS220"};
        for (String genomeId : genomeIds) {
            String result = exec.execute("genome " + genomeId);
            assertEquals("OK", result);
            assertEquals(genomeId, GenomeManager.getInstance().getCurrentGenome().getId());
        }
    }

    @Test
    public void testLoadGenomeFile() throws Exception {
        String[] genomePaths = new String[]{TestUtils.DATA_DIR + "genomes/hg18.unittest.genome"};
        String[] genomeIds = new String[]{"hg18.unittest"};
        int ind = 0;
        for (String genomePath : genomePaths) {
            String result = exec.execute("genome " + genomePath);
            assertEquals("OK", result);
            assertEquals(genomeIds[ind++], GenomeManager.getInstance().getCurrentGenome().getId());
        }
    }

    @Test
    public void testLoadGenomeFastaFile() throws Exception {
        String[] genomePaths = new String[]{TestUtils.DATA_DIR + "fasta/ecoli_out.padded.fasta"};
        String[] genomeIds = genomePaths;
        int ind = 0;
        for (String genomePath : genomePaths) {
            String result = exec.execute("genome " + genomePath);
            assertEquals("OK", result);
            assertEquals(genomeIds[ind++], GenomeManager.getInstance().getCurrentGenome().getId());
        }
    }

    @Test
    public void testLoadGenomesFail() throws Exception {
        String startId = genome.getId();
        String[] genomeIds = new String[]{"hg1920", "inch10", "doctor5"};
        for (String genomeId : genomeIds) {
            String result = exec.execute("genome " + genomeId);
            assertTrue(result.toLowerCase().startsWith(("error")));
            assertEquals(startId, GenomeManager.getInstance().getCurrentGenome().getId());
        }
    }

}
