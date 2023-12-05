package org.broad.igv.ucsc.bb;

import org.broad.igv.data.DataTile;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.IOException;
import java.util.List;

import static org.junit.Assert.*;

public class BBDataSourceTest {

    @Test
    public void testUncompressedBigwig() throws IOException {

        //chr21:19,146,376-19,193,466
        String url = "https://s3.amazonaws.com/igv.org.test/data/uncompressed.bw";
        String chr = "chr21";
        int start = 0;
        int end = Integer.MAX_VALUE;
        double bpPerPixel = 6191354.824;

        BBFile bbFile = new BBFile(url, null);

        assertNotNull(bbFile);

//        const bwReader = new BWReader({url: url})
//        const features = await bwReader.readFeatures(chr, start, chr, end, bpPerPixel)
//        assert.equal(features.length, 8)   // Verified in iPad app

    }

    @Test
    public void testBigWigZoom() throws IOException {

        //chr21:19,146,376-19,193,466
        String url = "https://www.encodeproject.org/files/ENCFF000ARZ/@@download/ENCFF000ARZ.bigWig";

        //String url = TestUtils.DATA_DIR + "bb/fixedStep.bw";
        String chr = "1";

        Genome genome = TestUtils.mockUCSCGenome();
        BBFile bbFile = new BBFile(url, genome);
        BBDataSource bbDataSource = new BBDataSource(bbFile, genome);

        int zoomLevel = 1;
        List<LocusScore> scores = bbDataSource.getPrecomputedSummaryScores(chr, 0, Integer.MAX_VALUE, zoomLevel);
        assertNotNull(scores);
        assertTrue(scores.size() > 0);

        //High resolutions -- there should be no precomputed scores (i.e. no zoom data).

        zoomLevel = 20;
        scores = bbDataSource.getPrecomputedSummaryScores(chr, 0, Integer.MAX_VALUE, zoomLevel);
        assertNull(scores);

    }

    @Test
    public void testBigWigWig() throws IOException {

        //chr21:19,146,376-19,193,466
        String url = "https://www.encodeproject.org/files/ENCFF000ARZ/@@download/ENCFF000ARZ.bigWig";

        //String url = TestUtils.DATA_DIR + "bb/fixedStep.bw";
        String chr = "1";
        int start = 72464570;
        int end = 72464687;

        // Expected values -- from manual inspection at igv.org/app
        float [] expectedValues =  {1.0f, 1.96f, 3.0f, 3.0f, 4.16f, 6.0f};

        Genome genome = TestUtils.mockUCSCGenome();
        BBFile bbFile = new BBFile(url, genome);
        BBDataSource bbDataSource = new BBDataSource(bbFile, genome);

        DataTile data = bbDataSource.getRawData(chr, start, end);
        assertNotNull(data);
        assertTrue(data.getValues().length == expectedValues.length);
        for(int i=0; i<expectedValues.length; i++) {
            assertEquals(expectedValues[i], data.getValues()[i], 0.000000001f);
        }

    }
}