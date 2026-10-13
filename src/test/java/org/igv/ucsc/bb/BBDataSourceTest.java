package org.igv.ucsc.bb;

import org.igv.data.DataTile;
import org.igv.feature.LocusScore;
import org.igv.feature.genome.Genome;
import org.igv.util.TestUtils;
import org.junit.Test;

import java.io.IOException;
import java.util.List;

import static org.junit.Assert.*;

public class BBDataSourceTest {


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

    /**
     * Test a big-endian bigWig.  Contents:
     * chr1  bedGraph section   100-200 1.5, 300-400 2.5
     * chr2  fixedStep section  1000-1050 0.25, 1100-1150 0.5, 1200-1250 0.75
     * One zoom level, reduction 1000.
     */
    @Test
    public void testBigEndian() throws IOException {

        String path = TestUtils.DATA_DIR + "bb/bigEndian.bw";

        Genome genome = TestUtils.mockUCSCGenome();
        BBFile bbFile = new BBFile(path, genome);
        assertEquals(2, bbFile.getChromosomeCount());
        BBDataSource bbDataSource = new BBDataSource(bbFile, genome);

        DataTile data = bbDataSource.getRawData("chr1", 0, 1000);
        assertArrayEquals(new int[]{100, 300}, data.getStartLocations());
        assertArrayEquals(new int[]{200, 400}, data.getEndLocations());
        assertArrayEquals(new float[]{1.5f, 2.5f}, data.getValues(), 0);

        // chr2 has chromosome id 1, which reads as 16777216 in the wrong byte order
        data = bbDataSource.getRawData("chr2", 0, 2000);
        assertArrayEquals(new int[]{1000, 1100, 1200}, data.getStartLocations());
        assertArrayEquals(new int[]{1050, 1150, 1250}, data.getEndLocations());
        assertArrayEquals(new float[]{0.25f, 0.5f, 0.75f}, data.getValues(), 0);

        List<LocusScore> scores = bbDataSource.getPrecomputedSummaryScores("chr2", 0, Integer.MAX_VALUE, 0);
        assertEquals(1, scores.size());
        assertEquals(1000, scores.get(0).getStart());
        assertEquals(1250, scores.get(0).getEnd());
        assertEquals(0.5f, scores.get(0).getScore(), 0);
    }
}