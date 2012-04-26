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

import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.tools.IgvTools;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.tribble.Feature;
import org.broad.tribble.TribbleException;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertNotNull;


/**
 * @author Jim Robinson
 * @date 10/3/11
 */
public class TrackLoaderTest {

    TrackLoader trackLoader;

    @Before
    public void setUp() throws Exception {
        TestUtils.setUpHeadless();
        trackLoader = new TrackLoader();
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
        TestUtils.stopGUI();
    }

    @Test
    public void testLoadBEDIndexed() throws Exception {
        String filepath = TestUtils.DATA_DIR + "bed/intervalTest.bed";
        TestUtils.createIndex(filepath);
        tstLoadFi(filepath, 1);

    }

    @Test
    public void testLoadBEDNotIndexed() throws Exception {
        String filepath = TestUtils.DATA_DIR + "bed/intervalTest.bed";
        if (TrackLoader.isIndexed(filepath, null)) {
            File f = new File(filepath + ".idx");
            f.delete();
        }
        tstLoadFi(filepath, 1);

    }


    @Test(expected = TribbleException.MalformedFeatureFile.class)
    public void testBEDCodec1() throws Exception {
        String filepath = TestUtils.DATA_DIR + "bed/NA12878.deletions.10kbp.het.gq99.hand_curated.hg19.bed";
        tstLoadFi(filepath, null);
    }

    @Test
    public void testLoadSIF() throws Exception {
        String filepath = TestUtils.DATA_DIR + "sample/BRCA_sif.txt";
        //Sample information file, shouldn't have tracks. Not a great test
        tstLoadFi(filepath, 0);
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
        IGVFeature feat0 = (IGVFeature) features.get(0);
        IGVFeature feat1 = (IGVFeature) features.get(1);

        assertEquals(707 - 1, feat0.getStart());
        assertEquals(943, feat0.getEnd());
        assertEquals(7563 - 1, feat1.getStart());
        assertEquals(7938, feat1.getEnd());
    }

    @Test
    public void testLoadGFFAliasedChrs() throws Exception{
        String filepath = TestUtils.DATA_DIR + "gff/aliased.gff";
        TrackLoader loader = new TrackLoader();
        Genome genome = IgvTools.loadGenome(TestUtils.DATA_DIR + "genomes/hg18_truncated_aliased.genome", true);
        List<Track> tracks = loader.load(new ResourceLocator(filepath), genome);
        assertEquals(1, tracks.size());
        FeatureTrack track = (FeatureTrack) tracks.get(0);

        assertEquals("aliased.gff", track.getName());

        List<Feature> features = track.getFeatures("chr1", 0, Integer.MAX_VALUE);
        assertEquals(56, features.size());

        features = track.getFeatures("chr5", 0, Integer.MAX_VALUE);
        assertEquals(16, features.size());

        //Non-aliased
        features = track.getFeatures("NC_007072.3", 0, Integer.MAX_VALUE);
        assertEquals(30, features.size());
    }





    private List<Track> tstLoadFi(String filepath, Integer expected_tracks) throws Exception {
        Genome genome = TestUtils.loadGenome();
        return tstLoadFi(filepath, expected_tracks, genome);
    }

    private List<Track> tstLoadFi(String filepath, Integer expected_tracks, Genome genome) throws Exception {
        ResourceLocator locator = new ResourceLocator(filepath);

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
        tstLoadFi(TestUtils.DATA_DIR + "bed/canFam2_alias.bed", 1);
        String[] aliases = new String[]{"AAAA", "BBB", "CCC"};
        for (String alias : aliases) {
            Feature feat = FeatureDB.getFeature(alias);
            assertNotNull(feat);
        }

    }

    private static String[] filenames = new String[]{"/bb/chr21.refseq.bb", "/bed/MT_test.bed", "/bed/Unigene.sample.bed",
            "/bed/test.bed", "/cn/HindForGISTIC.hg16.cn", "/folder with spaces/test.wig",
            "/gct/igv_test2.gct", "/gct/affy_human_mod.gct", "/gff/gene.gff3", "/igv/MIP_44.cn",//"/gc/chr1.txt",
            "/maf/TCGA_GBM_Level3_Somatic_Mutations_08.28.2008.maf.gz", "/psl/fishBlat.psl", "/sam/test_2.sam",
            "/seg/canFam2_hg18.seg", "/wig/test.wig"};

    @Test
    public void testFilesHeadless() throws Exception {
        Genome genome = TestUtils.loadGenome();
        for (String finame : filenames) {
            tstLoadFi(TestUtils.DATA_DIR + finame, null, genome);
        }
    }

    @Test
    public void testFilesHeaded() throws Exception {
        TestUtils.startGUI();

        String ex_filename = "/vcf/example4-last-gsnap-2_fixed.vcf";
        Genome genome = TestUtils.loadGenome();
        List<String> finames = new ArrayList<String>(Arrays.asList(filenames));

        finames.add(ex_filename);

        for (String finame : finames) {
            tstLoadFi(TestUtils.DATA_DIR + finame, null, genome);
        }
        TestUtils.stopGUI();
    }


}
