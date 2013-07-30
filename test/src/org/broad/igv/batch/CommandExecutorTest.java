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

import org.apache.commons.lang.StringUtils;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.dev.api.batch.Command;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.RegionScoreType;
import org.broad.igv.track.Track;
import org.broad.igv.ui.AbstractHeadedTest;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.IGVTestHeadless;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.*;
import java.util.List;

import static junit.framework.Assert.*;

/**
 * User: jacob
 * Date: 2012/03/21
 */
public class CommandExecutorTest extends AbstractHeadedTest {

    CommandExecutor exec = new CommandExecutor();
    private final String snapshotDir = TestUtils.TMP_OUTPUT_DIR;

    @Rule
    public TestRule testTimeout = new Timeout((int) 1800000);

    @Before
    public void setUp() throws Exception {
        super.setUp();
        Globals.setBatch(true);
        igv.loadGenome(TestUtils.defaultGenome, null);
        igv.newSession();
        exec.setSnapshotDirectory(snapshotDir);
    }

    @Test
    public void testRegionNoname() throws Exception{
        tstRegion(null);
    }

    @Test
    public void testRegionName() throws Exception{
        tstRegion("myregion");
    }

    public void tstRegion(String desc) throws Exception {

        String regionStr = "chr1:50-1000";
        String descstr = desc != null ? desc : "";
        exec.execute("region " + regionStr + " " + descstr);

        Collection<RegionOfInterest> regions = IGV.getInstance().getSession().getAllRegionsOfInterest();
        assertEquals(1, regions.size());

        ArrayList<RegionOfInterest> regionsAL = (ArrayList<RegionOfInterest>) regions;
        RegionOfInterest region = regionsAL.get(0);
        assertEquals("chr1", region.getChr());
        assertEquals(49, region.getStart());
        assertEquals(1000, region.getEnd());
        assertEquals(descstr, region.getDescription());

    }

    @Test
    @Ignore
    public void stressTestSnapshotsHG00171() throws Exception{
        PreferenceManager.getInstance().put(PreferenceManager.SAM_MAX_VISIBLE_RANGE, "1000");

        String interv0 = "chr1:151666000-152666000";
        String interv1 = "chr1:154666000-155666000";
        String[] intervals = {interv0, interv1};

        String filePath = TestUtils.LARGE_DATA_DIR + "HG00171.hg18.bam";

        stressTstSnapshots(filePath, intervals);
    }

    @Test
    @Ignore
    public void stressTestSnapshotsBodymap() throws Exception{
        PreferenceManager.getInstance().put(PreferenceManager.SAM_DOWNSAMPLE_READS, "true");
        PreferenceManager.getInstance().put(PreferenceManager.SAM_SAMPLING_COUNT, "100");
        PreferenceManager.getInstance().put(PreferenceManager.SAM_MAX_VISIBLE_RANGE, "1000");

        String interv0 = "chr12:97,509,534-97,521,909"; //SLC25A3
        String interv1 = "chrX:153,366,844-153,374,196"; //SLC10A3
        String[] intervals = {interv0, interv1};

        String filePath = TestUtils.DATA_DIR + "sessions/bodymap_3tissue.xml";

        stressTstSnapshots(filePath, intervals);
    }



    /**
     * Take a large number of snapshots, make sure they all
     * actually show data.
     * @throws Exception
     */
    public void stressTstSnapshots(String filePath, String[] intervals) throws Exception{

        exec.execute("load " + filePath);
        //exec.execute(("setSleepInterval 10000"));

        //For each interval we base our expected size on the first snapshot
        Map<String, Long> intervalSizeMap = new HashMap<String, Long>(intervals.length);

        Long expSize;
        long margin;
        int numTrials = 50;
        for(int tri=0; tri < numTrials; tri++){

            int intInd = tri % intervals.length;
            String interval = intervals[intInd];
            expSize = intervalSizeMap.get(interval);

            exec.execute("goto " + interval);

            String outFileName = outFileBase + tri + ".png";
            File outFile = new File(snapshotDir, outFileName);
            if(outFile.exists()) outFile.delete();
            tstSnapshot(outFileName);

            long size = outFile.length();
            if(expSize == null){
                expSize = size;
                intervalSizeMap.put(interval, expSize);
            }
            margin = expSize / 10;
            long sizeDiff = Math.abs(size - expSize);
            //break;
            assertTrue(String.format("File size much different than expected. Trial %d, Diff = %d, margin = %d", tri, sizeDiff, margin),  sizeDiff < margin);
        }

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

    public static final String urlPathSpaces = "ftp://ftp.broadinstitute.org/distribution/igv/TEST/cpgIslands%20with%20spaces.hg18.bed";

    public static final String dirPathSpaces = TestUtils.DATA_DIR + "folder with spaces";
    public static final String fileName01 = "test.wig";
    public static final String fileNamePerc = "%test%2D.wig";
    public static final String fileNamePlus = "test+wp.wig";


    @Test
    public void testLoadFileSpaces() throws Exception{
        tstLoadFileSpaces(fileName01);
    }

    @Test
    public void testLoadFileSpacesPerc() throws Exception{
        tstLoadFileSpaces(fileNamePerc);
    }

    @Test
    public void testLoadFileSpacesPlus() throws Exception{
        tstLoadFileSpaces(fileNamePlus);
    }


    @Test
    public void testLoadFileURLSpaces() throws Exception{
        tstLoadFileURLSpaces(fileName01);
    }

    @Test
    public void testLoadFileURLSpacesPerc() throws Exception{
        tstLoadFileURLSpaces(fileNamePerc);
    }

    @Test
    public void testLoadFileURLSpacesPlus() throws Exception{
        tstLoadFileURLSpaces(fileNamePlus);
    }

    private void tstLoadFileURLSpaces(String filename) throws Exception{
        String fileURL = "file://" + org.broad.igv.util.StringUtils.encodeURL(new File(dirPathSpaces, filename).getAbsolutePath());
        exec.execute("load " + fileURL);
        TestUtils.assertTrackLoaded(IGV.getInstance(), filename);
    }

    private void tstLoadFileSpaces(String filename) throws Exception{
        File file = new File(dirPathSpaces, filename);
        exec.execute("load \"" + file.getPath() + "\"");
        TestUtils.assertTrackLoaded(IGV.getInstance(), filename);
    }

    @Test
    public void testLoadURL() throws Exception {
        //This is mostly to ruggedize test setup. The setup may load
        //reference/sequence tracks, we'd like to be able to change
        //test setup and not worry about this test.
        int beginTracks = igv.getAllTracks().size();

        String urlPath = urlPathSpaces;
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

    @Test
    public void testPreference() throws Exception{
        String key = PreferenceManager.DATA_SERVER_URL_KEY;
        String val = "myDataServerURL";

        assertNotSame(val, PreferenceManager.getInstance().getDataServerURL());

        exec.execute(String.format("preference %s %s", key, val));

        assertEquals(val, PreferenceManager.getInstance().getDataServerURL());
    }

    @Test
    public void testSnapshotsize() throws Exception{
        String filePath = TestUtils.DATA_DIR + "bam/NA12878.SLX.sample.bam";
        int numLoads = 1;

        for(int ii = 0; ii < numLoads; ii++){
            IGV.getInstance().loadResources(Arrays.asList(new ResourceLocator(filePath)));
        }
        exec.execute("goto chr1:9,713,386-9,733,865");


        int minHeight = (int) Toolkit.getDefaultToolkit().getScreenSize().getHeight() - 150;
        int maxHeight = minHeight + 200;

        String outFileName =  minHeight + ".png";

        exec.execute("maxpanelheight " + maxHeight);
        exec.execute("snapshot " + outFileName);

        File outputFile = new File(snapshotDir, outFileName);
        BufferedImage image = ImageIO.read(outputFile);

        assertTrue("Output image height " + image.getHeight() + " is not at least " + minHeight, image.getHeight() > minHeight);

        int remAlphaMask = 0x00ffffff;

        int numBlackPix = 0;
        for(int yy = image.getMinY(); yy < image.getHeight(); yy++){
            for(int xx = image.getMinX(); xx < image.getWidth(); xx++){
                int color = image.getRGB(xx, yy) & remAlphaMask;
                numBlackPix += color == 0 ? 1 : 0;
            }
        }

        //Just making sure we don't trivially satisfy the problem
        assertTrue(numBlackPix > 100);

        int totalPix = image.getHeight()*image.getWidth();
        assertTrue("Too much of the snapshot is black", numBlackPix < totalPix * 0.1);
    }

    @Test
    public void testCustomCommand_Echo() throws Exception {
        String cmd = EchoCommand.class.getName();
        String otherArgs = "fly high free bird";
        String fullCmd = String.format("%s %s", cmd, otherArgs);
        String response = exec.execute(fullCmd);

        assertEquals(otherArgs, response);
    }

    public static class EchoCommand implements Command {
        @Override
        public String run(List<String> args) {
            return StringUtils.join(args, " ");
        }
    }
}
