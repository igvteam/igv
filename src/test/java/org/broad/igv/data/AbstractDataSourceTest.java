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
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackType;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Dec 19, 2009
 * Time: 7:17:32 PM
 * To change this template use File | Settings | File Templates.
 */
public class AbstractDataSourceTest extends AbstractHeadlessTest {

    /**
     * # When viewed as a heatmap this feature shows a 1-pixel break at position 23,314,405  when the
     * # zoom is at  chr20:23,314,141-23,314,667
     * chr20   23313029    23316433    91.1
     */
    @Test
    public void testRT_134467() {
        int[] starts = {23313029};
        int[] ends = {23316433};
        float[] values = {91.1f};

        int s = 23313929;
        int e = 23314405;

        TestDataSource ds = new TestDataSource(starts, ends, values);

        SummaryTile tile = ds.computeSummaryTile("chr20", s, e, 1);
        List<LocusScore> scores = tile.getScores();

        // Scores should be within 10 +/- 0.5,  and the mean should be very close to 10

        for (LocusScore score : scores) {
            float v = score.getScore();
            assertEquals(91.1f, v, 0.00001);
        }

        assertEquals(1, scores.size());


    }

    @Test
    public void testGetSummaryScoresForRange() {

        TestDataSource ds = new TestDataSource();

        SummaryTile tile = ds.computeSummaryTile("", 0, 10000, 1);

        List<LocusScore> scores = tile.getScores();

        // Scores should be within 10 +/- 0.5,  and the mean should be very close to 10
        float sum = 0.0f;
        long totPoints = 0;
        for (LocusScore score : scores) {
            float v = score.getScore();
            assertTrue((v >= 9.5f && v <= 10.5f));
            int numPoints = score.getEnd() - score.getStart();
            sum += numPoints * v;
            totPoints += numPoints;
        }
        double mean = sum / totPoints;
        assertEquals(10.0, mean, 1.0e-2);
    }

    @Test
    public void testGetSummaryScoresForSNPs() throws Exception {

        ResourceLocator locator = new ResourceLocator(TestUtils.DATA_DIR + "cn/multi_snp.cn");
        Genome genome = TestUtils.loadGenome();
        IGVDataset ds = new IGVDataset(locator, genome);
        DatasetDataSource dataSource = new DatasetDataSource("Sample1", ds, genome);
        String chr = "chr10";

        dataSource.cacheSummaryTiles = false;
        int zreq = 22;
        int half_width = 20;
        int[] starts = {72644150, 72698871, 72729621, 89614266, 89614367, 89614406, 89614478};
        for (int start : starts) {
            int end = start + 1 + half_width;
            start -= half_width;
            List<LocusScore> scores = dataSource.getSummaryScoresForRange(chr, start, end, zreq);
            assertEquals(1, scores.size());
        }
        int start = starts[0] - 100;
        int end = starts[starts.length - 1] + 100;
        List<LocusScore> scores = dataSource.getSummaryScoresForRange(chr, start, end, 22);
        //The last few get combined into 1 tile
        assertEquals(4, scores.size());

    }


    public class TestDataSource extends AbstractDataSource {


        private int nPts;
        int[] starts;
        int[] ends;
        float[] values;
        String[] probes;

        TestDataSource(int[] starts, int[] ends, float[] values) {
            super(null);
            nPts = starts.length;
            this.starts = starts;
            this.ends = ends;
            this.values = values;
            this.probes = null;
        }

        TestDataSource() {
            super(null);
            nPts = 10000;
            starts = new int[nPts];
            ends = new int[nPts];
            values = new float[nPts];
            probes = new String[nPts];
            for (int i = 0; i < nPts; i++) {
                starts[i] = i;
                ends[i] = i + 1;
                values[i] = (float) (9.5 + Math.random());
                probes[i] = "probe_" + i;
            }


        }

        protected DataTile getRawData(String chr, int startLocation, int endLocation) {
            return new DataTile(starts, ends, values, probes);
        }

        @Override
        protected List<LocusScore> getPrecomputedSummaryScores(String chr, int startLocation, int endLocation, int zoom) {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }


        @Override
        public int getLongestFeature(String chr) {
            return 1000;
        }

        public double getMedian(int zoom, String chr) {
            return 0;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public double getDataMax() {
            return 0;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public double getDataMin() {
            return 0;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public TrackType getTrackType() {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }
    }
}
