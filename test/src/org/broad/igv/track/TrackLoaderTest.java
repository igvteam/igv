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

package org.broad.igv.track;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.tools.IgvTools;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.igv.variant.VariantTrack;
import htsjdk.tribble.Feature;
import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

import java.io.File;
import java.util.*;

import static junit.framework.Assert.*;


/**
 * @author Jim Robinson
 * @date 10/3/11
 */
public class TrackLoaderTest extends AbstractHeadlessTest {

    TrackLoader trackLoader;

    @Rule
    public TestRule testTimeout = new Timeout((int) (500 * 60e3));

    @Before
    public void setUp() throws Exception {
        super.setUp();
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
    public void testLoadBEDFtp() throws Exception {
        String filepath = "ftp://ftp.broadinstitute.org/distribution/igv/TEST/cpgIslands with spaces.hg18.bed";
        tstLoadFi(filepath, 1, false);
    }

    @Test
    public void testLoadBEDNotIndexed() throws Exception {
        String filepath = TestUtils.DATA_DIR + "bed/intervalTest.bed";
        if (TrackLoader.isIndexed(new ResourceLocator(filepath), null)) {
            File f = new File(filepath + ".idx");
            f.delete();
        }
        tstLoadFi(filepath, 1, false);

    }


    @Test
    public void testBEDCodec1() throws Exception {
        // This file has an extra tab in the second to last line.  Should have no effect, but earlier versions
        // of IGV would split on tabs only and throw a runtime exception.
        String filepath = TestUtils.DATA_DIR + "bed/multi_tab.bed";
        tstLoadFi(filepath, null, false);
    }

    /**
     * Make sure we throw an exception when parsing a bad file
     * @throws Exception
     */
    @Test(expected = RuntimeException.class)
    public void testBEDCodec2() throws Exception {
        String filepath = TestUtils.DATA_DIR + "bed/bad_file.bed";
        tstLoadFi(filepath, null, false);
    }

    @Test
    public void testLoadSIF() throws Exception {
        String filepath = TestUtils.DATA_DIR + "sample/BRCA_sif.txt";
        List<String> expectedAttributeNames = Arrays.asList("TCGA_EXPERIMENT", "TCGA_BATCH", "TUMOR_NORMAL",
                "BIRDSEED_GENDER", "LEVEL2_NOISE", "LEVEL3_SEGMENT_COUNT", "PURITY	PLOIDY", "DELTA",
                "CANCER_DNA_FRACTION", "SUBCLONAL_GENOME_FRACTION");

        tstLoadFi(filepath, 0, false);  //Sample information file, shouldn't have tracks.

        Set<String> attrNames = new HashSet<String>(AttributeManager.getInstance().getAttributeNames());

        assertTrue(attrNames.size() >= expectedAttributeNames.size());  // Can be larger because of default attributes
        for (String name : expectedAttributeNames) {
            assertTrue(expectedAttributeNames.contains(name));
        }

    }

    @Test
    public void testLoadGFF() throws Exception {
        String filepath = TestUtils.DATA_DIR + "gff/simfeatures.gff3";
        List<Track> tracks = trackLoader.load(new ResourceLocator(filepath), genome);
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
    public void testLoadGFFAliasedChrs() throws Exception {
        String filepath = TestUtils.DATA_DIR + "gff/aliased.unsorted.gff";
        Genome genome = IgvTools.loadGenome(TestUtils.DATA_DIR + "genomes/hg18_truncated_aliased.genome");
        List<Track> tracks = trackLoader.load(new ResourceLocator(filepath), genome);
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

    /**
     * Test that we properly load a vcf file, even though the
     * extension is upper-case. IGV-2012
     *
     * @throws Exception
     */
    @Test
    public void testLoadVCFUpperCase() throws Exception {
        String filepath = TestUtils.DATA_DIR + "vcf/HC_MOD_CAPS.VCF";
        ResourceLocator locator = new ResourceLocator(filepath);
        TestUtils.createIndex(filepath);

        List<Track> tracks = trackLoader.load(locator, genome);
        assertEquals(1, tracks.size());
        assertTrue("VCF file loaded incorrect track type", tracks.get(0) instanceof VariantTrack);

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

    static String[] filenamesNoIndex = new String[]{"cufflinks/sample_genes.fpkm_tracking", "cufflinks/sample_gene_exp.diff", "/bb/chr21.refseq.bb", "/bed/MT_test.bed", "/bed/Unigene.sample.bed",
            "/bed/test.bed", "/cn/HindForGISTIC.hg16.cn", "/folder with spaces/test.wig",
            "/gct/igv_test2.gct", "/gct/affy_human_mod.gct", "/gff/gene.unsorted.gff3", "/igv/MIP_44.cn",
            "/maf/TCGA_GBM_Level3_Somatic_Mutations_08.28.2008.maf.gz",
            "/psl/fishBlat.psl", "/seg/canFam2_hg18.seg", "/wig/test.wig"};

    private static ArrayList<String> tmp = new ArrayList<String>(filenamesNoIndex.length);
    static String[] filenamesTryIndex;

    static {
        tmp.add("/sam/test_2.sam");
        filenamesTryIndex = tmp.toArray(new String[tmp.size()]);
    }

    public void tstFilesHeadless(String[] filenames, boolean tryIndex) throws Exception {
        Genome genome = TestUtils.loadGenome();
        for (String finame : filenames) {
            tstLoadFi(TestUtils.DATA_DIR + finame, null, genome, tryIndex);
        }
    }

    @Test
    public void testFilesHeadlessTryIndex() throws Exception {
        tstFilesHeadless(filenamesTryIndex, true);
    }

    @Test
    public void testFilesHeadlessNoIndex() throws Exception {
        tstFilesHeadless(filenamesNoIndex, false);
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

    @Test
    public void testLoadBAMFtp() throws Exception {
        String path = "ftp://ftp.broadinstitute.org/pub/igv/TEST/HG00171.hg18.bam";
        tstLoadFi(path, 2, false);
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
