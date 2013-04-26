/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.data;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.data.seg.SegmentFileParser;
import org.broad.igv.data.seg.SegmentedAsciiDataSet;
import org.broad.igv.data.seg.SegmentedDataSet;
import org.broad.igv.data.seg.SegmentedDataSource;
import org.broad.igv.feature.LocusScore;
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
        float[] expScores = new float[]{2,  2,  3,  3,  4,  4,  5,  5,  1,  1};
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

        return new CombinedDataSource(sourceA, sourceB, operation);
    }

    private SegmentedAsciiDataSet getSegDataSet(String path){
        ResourceLocator locator = new ResourceLocator(path);
        SegmentFileParser parser = new SegmentFileParser(locator);
        return parser.loadSegments(locator, genome);
    }
}
