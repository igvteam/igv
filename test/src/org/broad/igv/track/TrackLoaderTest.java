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

package org.broad.igv.track;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.tools.IgvTools;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.Feature;
import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

import java.io.File;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static junit.framework.Assert.*;


/**
 * @author Jim Robinson
 * @date 10/3/11
 */
public class TrackLoaderTest extends AbstractHeadlessTest {

    TrackLoader trackLoader;

    @Rule
    public TestRule testTimeout = new Timeout((int) (5 * 60e3));

    @Before
    public void setUp() throws Exception {
        trackLoader = new TrackLoader();
    }

    @Test
    public void testLoadBEDIndexed() throws Exception {
        String filepath = TestUtils.DATA_DIR + "bed/intervalTest.bed";
        TestUtils.createIndex(filepath);
        tstLoadFi(filepath, 1, true);
    }

    @Test
    public void testLoadBEDtxt() throws Exception {
        String filepath = TestUtils.DATA_DIR + "bed/intervalTest.bed.txt";
        TestUtils.createIndex(filepath);
        tstLoadFi(filepath, 1, true);
    }

    @Test
    public void testLoadBEDNotIndexed() throws Exception {
        String filepath = TestUtils.DATA_DIR + "bed/intervalTest.bed";
        if (TrackLoader.isIndexed(filepath, null)) {
            File f = new File(filepath + ".idx");
            f.delete();
        }
        tstLoadFi(filepath, 1, false);

    }


    @Test(expected = RuntimeException.class)
    public void testBEDCodec1() throws Exception {
        String filepath = TestUtils.DATA_DIR + "bed/NA12878.deletions.10kbp.het.gq99.hand_curated.hg19.bed";
        tstLoadFi(filepath, null, false);
    }

    @Test
    public void testLoadSIF() throws Exception {
        String filepath = TestUtils.DATA_DIR + "sample/BRCA_sif.txt";
        //Sample information file, shouldn't have tracks. Not a great test
        tstLoadFi(filepath, 0, false);
    }

    @Test
    public void testLoadGFF() throws Exception {
        String filepath = TestUtils.DATA_DIR + "gff/simfeatures.gff3";
        TrackLoader loader = new TrackLoader();
        List<Track> tracks = loader.load(new ResourceLocator(filepath), TestUtils.loadGenome());
        assertEquals(1, tracks.size());
        FeatureTrack track = (FeatureTrack) tracks.get(0);

        assertEquals("notmeaningful", track.getName());

        List<Feature> features = track.getFeatures("chr1", 0, Integer.MAX_VALUE);
        assertEquals(2, features.size());
        //Since the features are type exon, they will be added as exons
        IGVFeature feat0 = ((IGVFeature) features.get(0)).getExons().get(0);
        IGVFeature feat1 = ((IGVFeature) features.get(1)).getExons().get(0);

        assertEquals(707 - 1, feat0.getStart());
        assertEquals(943, feat0.getEnd());
        assertEquals(7563 - 1, feat1.getStart());
        assertEquals(7938, feat1.getEnd());
    }

    @Test
    public void testLoadGFFAliasedChrs() throws Exception {
        String filepath = TestUtils.DATA_DIR + "gff/aliased.unsorted.gff";
        TrackLoader loader = new TrackLoader();
        Genome genome = IgvTools.loadGenome(TestUtils.DATA_DIR + "genomes/hg18_truncated_aliased.genome", true);
        List<Track> tracks = loader.load(new ResourceLocator(filepath), genome);
        assertEquals(1, tracks.size());
        FeatureTrack track = (FeatureTrack) tracks.get(0);

        assertEquals("aliased.unsorted.gff", track.getName());

        List<Feature> features = track.getFeatures("chr1", 0, Integer.MAX_VALUE);
        assertEquals(56, features.size());

        features = track.getFeatures("chr5", 0, Integer.MAX_VALUE);
        assertEquals(16, features.size());

        //Non-aliased
        features = track.getFeatures("NC_007072.3", 0, Integer.MAX_VALUE);
        assertEquals(30, features.size());
    }


    private List<Track> tstLoadFi(String filepath, Integer expected_tracks, boolean makeIndex) throws Exception {
        Genome genome = TestUtils.loadGenome();
        return tstLoadFi(filepath, expected_tracks, genome, makeIndex);
    }

    private List<Track> tstLoadFi(String filepath, Integer expected_tracks, Genome genome, boolean makeIndex) throws Exception {
        return tstLoadFi(this.trackLoader, filepath, expected_tracks, genome, makeIndex);
    }

    static List<Track> tstLoadFi(TrackLoader trackLoader, String filepath, Integer expected_tracks, Genome genome, boolean makeIndex) throws Exception {
        ResourceLocator locator = new ResourceLocator(filepath);

        //Try creating an index
        //UI would ask for confirmation
        if (makeIndex) {
            try {
                TestUtils.createIndex(filepath);
            } catch (Exception e) {

            }
        }
        List<Track> tracks = trackLoader.load(locator, genome);
        if (expected_tracks != null) {
            assertEquals(expected_tracks.intValue(), tracks.size());
            if (expected_tracks == 0) {
                return tracks;
            }
        }

        Track track = tracks.get(0);
        assertEquals(locator, track.getResourceLocator());

        return tracks;
    }

    @Test
    public void testBEDLoadsAliases() throws Exception {
        tstLoadFi(TestUtils.DATA_DIR + "bed/canFam2_alias.bed", 1, false);
        String[] aliases = new String[]{"AAAA", "BBB", "CCC"};
        for (String alias : aliases) {
            Feature feat = FeatureDB.getFeature(alias);
            assertNotNull(feat);
        }

    }

    static String[] filenames = new String[]{"/bb/chr21.refseq.bb", "/bed/MT_test.bed", "/bed/Unigene.sample.bed",
            "/bed/test.bed", "/cn/HindForGISTIC.hg16.cn", "/folder with spaces/test.wig",
            "/gct/igv_test2.gct", "/gct/affy_human_mod.gct", "/gff/gene.unsorted.gff3", "/igv/MIP_44.cn",//"/gc/chr1.txt",
            "/maf/TCGA_GBM_Level3_Somatic_Mutations_08.28.2008.maf.gz", "/psl/fishBlat.psl", "/sam/test_2.sam",
            "/seg/canFam2_hg18.seg", "/wig/test.wig"};

    @Test
    public void testFilesHeadless() throws Exception {
        Genome genome = TestUtils.loadGenome();
        for (String finame : filenames) {
            tstLoadFi(TestUtils.DATA_DIR + finame, null, genome, true);
        }
    }

    @Test
    public void testWigAndBigWig() throws Exception {
        String wigPath = TestUtils.DATA_DIR + "wig/test_fixedStep.wig";
        String bigWigPath = TestUtils.DATA_DIR + "wig/test_fixedStep.bigwig";
        List<Track> wigtracks = tstLoadFi(wigPath, 1, false);
        List<Track> bigWigtracks = tstLoadFi(bigWigPath, 1, false);

        String chr = "chr19";
        int start = 5930725;
        int end = 5930764;
        int zoom = 21;

        DataSourceTrack wigTrack = (DataSourceTrack) wigtracks.get(0);
        DataSourceTrack bigWigTrack = (DataSourceTrack) bigWigtracks.get(0);

        int trials = 10;

        for (int ii = 0; ii < trials; ii++) {
            int strt = start + ii;
            List<LocusScore> wigScores = wigTrack.getSummaryScores(chr, strt, end, zoom);
            List<LocusScore> bigWigScores = bigWigTrack.getSummaryScores(chr, strt, end, zoom);
            assertEquals(wigScores.size(), bigWigScores.size());
            int ind = 0;
            for (LocusScore ws : wigScores) {
                LocusScore bws = bigWigScores.get(ind);
                assertEquals(ws.getScore(), bws.getScore());
                ind++;
            }
        }

    }

    @Test
    public void testBedAndBigBed() throws Exception {
        String bedPath = TestUtils.DATA_DIR + "bed/Unigene.sample.nolong.bed";
        String bigBedPath = TestUtils.DATA_DIR + "bed/Unigene.sample.nolong.bigbed";

        File bigBedFile = new File(bigBedPath);

        //Need to index so query of bed file is accurate
        TestUtils.createIndex(bedPath);

        List<Track> bedtracks = tstLoadFi(bedPath, 1, false);
        List<Track> bigBedtracks = tstLoadFi(bigBedPath, 1, false);


        String chr = "chr2";
        int start = 178711404 - 1;
        int end = 179189619 + 1;

        FeatureTrack bedTrack = (FeatureTrack) bedtracks.get(0);
        FeatureTrack bigBedTrack = (FeatureTrack) bigBedtracks.get(0);

        //Multiple trials because we're concerned about file open/close issues
        int trials = 10;

        for (int ii = 0; ii < trials; ii++) {
            int strt = start + ii;
            List<Feature> bedFeatures = bedTrack.getFeatures(chr, strt, end);
            List<Feature> bigBedFeatures = bigBedTrack.getFeatures(chr, strt, end);
            TestUtils.assertFeatureListsEqual(bedFeatures.iterator(), bigBedFeatures.iterator());

            //NOT FOOLPROOF
            assertTrue(bigBedFile.canWrite());
        }


    }

    /**
     * Test loading segmented data file from a sql database, using a profile
     *
     * @throws Exception
     */
    @Test
    public void testLoadSegProfile() throws Exception {
        String path = TestUtils.DATA_DIR + "sql/seg_canFam2_profile.dbxml";

        int expectedTracks = 6;
        List<Track> tracks = trackLoader.load(new ResourceLocator(path), genome);
        assertEquals(expectedTracks, tracks.size());
        Set<String> expSampleIds = new HashSet<String>(Arrays.asList("0123-A", "0123-B-1", "0123-C-1", "0123-C-2", "0123-C-3"));
        Set<String> actSampleIds = new HashSet<String>(5);
        for (Track track : tracks) {
            if (track instanceof DataSourceTrack) {
                actSampleIds.add(track.getName());
            }
        }
        assertEquals(expSampleIds, actSampleIds);
    }

    /**
     * Test loading sample information file from a sql database, using a profile
     *
     * @throws Exception
     */
    @Test
    public void testLoadSampleInfoProfile() throws Exception {

        AttributeManager.getInstance().clearAllAttributes();
        String path = TestUtils.DATA_DIR + "sql/sampleinfo_brca_sif_profile.dbxml";

        int expectedTracks = 0;
        List<Track> tracks = trackLoader.load(new ResourceLocator(path), genome);
        assertEquals(expectedTracks, tracks.size());

        String[] attrNames = "TCGA_EXPERIMENT	TCGA_BATCH	TUMOR_NORMAL	BIRDSEED_GENDER	LEVEL2_NOISE	LEVEL3_SEGMENT_COUNT	PURITY	PLOIDY	DELTA	CANCER_DNA_FRACTION	SUBCLONAL_GENOME_FRACTION".split("\\s+");
        Set<String> expAttrNames = new HashSet<String>(Arrays.asList(attrNames));
        List<String> actAttrNames = AttributeManager.getInstance().getAttributeNames();
        actAttrNames.remove("NAME");
        actAttrNames.remove("DATA TYPE");
        actAttrNames.remove("DATA FILE");
        assertEquals(actAttrNames.size(), expAttrNames.size());
        for (String attrName : actAttrNames) {
            assertTrue(expAttrNames.contains(attrName));
        }
    }

    //@Test
    public void testLoadScratch() throws Exception {
        String path = "http://www.broadinstitute.org/igvdata/annotations/hg19/EnsemblGenes.ensGene";
        tstLoadFi(path, 1, genome, false);
    }
}
