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

package org.broad.igv.data;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.data.seg.SegmentFileParser;
import org.broad.igv.data.seg.SegmentedAsciiDataSet;
import org.broad.igv.data.seg.SegmentedDataSet;
import org.broad.igv.data.seg.SegmentedDataSource;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.DataSourceTrack;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.util.List;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertNotNull;

/**
 * @author jacob
 * @date 2013-Apr-25
 */
public class CombinedDataSourceTest extends AbstractHeadlessTest {

    /**
     * Just test that nothing crashes
     * @throws Exception
     */
    @Test
    public void testNoExceptions() throws Exception{
        String chr = "chr1";
        int start = 0;
        int end = (int) 100e6;
        int zoom = 0;

        CombinedDataSource combinedSource = getDataSource(CombinedDataSource.Operation.ADD);
        List<LocusScore> combinedScores = combinedSource.getSummaryScoresForRange(chr, start, end, zoom);
        for(LocusScore score: combinedScores){
            assertNotNull(score);
        }
    }

    /**
     * Test combination when the two data sources have the same tiling boundaries
     * shouldn't need to split anything up, easy case
     * @throws Exception
     */
    @Test
    public void testSameBoundaries() throws Exception{
        String chr = "chr1";
        int start = 0;
        int end = 1000;
        int zoom = 0;

        CombinedDataSource combinedSource = getDataSource(CombinedDataSource.Operation.ADD);
        List<LocusScore> combinedScores = combinedSource.getSummaryScoresForRange(chr, start, end, zoom);

        int[] expStarts = new int[]{0,100,200,300};
        int[] expEnds = new int[]{100,200,300,400};
        float[] expScores = new float[]{1,2,3,4};
        int expCount = expStarts.length;

        assertEquals(expCount, combinedScores.size());

        int idx = 0;
        for(LocusScore score: combinedScores){
            assertEquals(expStarts[idx], score.getStart());
            assertEquals(expEnds[idx], score.getEnd());
            assertEquals(expScores[idx++], score.getScore());
        }
    }

    @Test
    public void testDifferentBoundaries() throws Exception{
        String chr = "chr2";
        int start = 0;
        int end = 1000;
        int zoom = 0;

        CombinedDataSource combinedSource = getDataSource(CombinedDataSource.Operation.ADD);
        List<LocusScore> combinedScores = combinedSource.getSummaryScoresForRange(chr, start, end, zoom);

        int[] expStarts =   new int[]{100,150,200,250,300,350,400,450,500,550};
        int[] expEnds =     new int[]{150,200,250,300,350,400,450,500,550,650};
        float[] expScores = new float[]{1,  2,  3,  3,  4,  4,  5,  5,  1,  1};
        int expCount = expStarts.length;

        assertEquals(expCount, combinedScores.size());

        int idx = 0;
        for(LocusScore score: combinedScores){
            assertEquals(expStarts[idx], score.getStart());
            assertEquals(expEnds[idx], score.getEnd());
            assertEquals(expScores[idx++], score.getScore());
        }
    }

    private CombinedDataSource getDataSource(CombinedDataSource.Operation operation){

        String pathA = TestUtils.DATA_DIR + "seg/toCombine_a.seg";
        String pathB = TestUtils.DATA_DIR + "seg/toCombine_b.seg";

        SegmentedDataSet dsA = getSegDataSet(pathA);
        SegmentedDataSet dsB = getSegDataSet(pathB);

        SegmentedDataSource sourceA = new SegmentedDataSource("0123-A", dsA);
        SegmentedDataSource sourceB = new SegmentedDataSource("0123-B-1", dsB);

        DataSourceTrack trackA = new DataSourceTrack(null, sourceA.getTrackIdentifier(), sourceA.getTrackIdentifier(), sourceA);
        DataSourceTrack trackB = new DataSourceTrack(null, sourceB.getTrackIdentifier(), sourceB.getTrackIdentifier(), sourceB);

        return new CombinedDataSource(trackA, trackB, operation);
    }

    private SegmentedAsciiDataSet getSegDataSet(String path){
        ResourceLocator locator = new ResourceLocator(path);
        SegmentFileParser parser = new SegmentFileParser(locator);
        return parser.loadSegments(locator, genome);
    }
}
