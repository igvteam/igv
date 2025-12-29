package org.igv.data;

import org.igv.AbstractHeadlessTest;
import org.igv.feature.LocusScore;
import org.igv.track.DataSourceTrack;
import org.igv.track.WindowFunction;
import org.igv.util.ResourceLocator;
import org.igv.util.TestUtils;
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
    public void testNoExceptions() throws Exception {
        String chr = "chr1";
        int start = 0;
        int end = (int) 100e6;
        int zoom = 0;

        CombinedDataSource combinedSource = getDataSource(CombinedDataSource.Operation.ADD);
        List<LocusScore> combinedScores = combinedSource.getSummaryScoresForRange(chr, start, end, zoom);
        for (LocusScore score : combinedScores) {
            assertNotNull(score);
        }
    }

    /**
     * Test combination when the two data sources have the same tiling boundaries
     * shouldn't need to split anything up, easy case
     * @throws Exception
     */
    @Test
    public void testSameBoundaries() throws Exception {
        String chr = "chr1";
        int start = 0;
        int end = 1000;
        int zoom = 0;

        CombinedDataSource combinedSource = getDataSource(CombinedDataSource.Operation.ADD);
        List<LocusScore> combinedScores = combinedSource.getSummaryScoresForRange(chr, start, end, zoom);

        int[] expStarts = new int[]{0, 100, 200, 300};
        int[] expEnds = new int[]{100, 200, 300, 400};
        float[] expScores = new float[]{1, 2, 3, 4};
        int expCount = expStarts.length;

        assertEquals(expCount, combinedScores.size());

        int idx = 0;
        for (LocusScore score : combinedScores) {
            assertEquals(expStarts[idx], score.getStart());
            assertEquals(expEnds[idx], score.getEnd());
            assertEquals(expScores[idx++], score.getScore());
        }
    }

    @Test
    public void testDifferentBoundaries() throws Exception {
        String chr = "chr2";
        int start = 0;
        int end = 1000;
        int zoom = 0;

        CombinedDataSource combinedSource = getDataSource(CombinedDataSource.Operation.ADD);
        List<LocusScore> combinedScores = combinedSource.getSummaryScoresForRange(chr, start, end, zoom);

        int[] expStarts = new int[]{100, 150, 200, 250, 300, 350, 400, 450, 500, 550};
        int[] expEnds = new int[]{150, 200, 250, 300, 350, 400, 450, 500, 550, 650};
        float[] expScores = new float[]{1, 2, 3, 3, 4, 4, 5, 5, 1, 1};
        int expCount = expStarts.length;

        assertEquals(expCount, combinedScores.size());

        int idx = 0;
        for (LocusScore score : combinedScores) {
            assertEquals(expStarts[idx], score.getStart());
            assertEquals(expEnds[idx], score.getEnd());
            assertEquals(expScores[idx++], score.getScore());
        }
    }

    private CombinedDataSource getDataSource(CombinedDataSource.Operation operation) {

        String pathA = TestUtils.DATA_DIR + "wig/toCombine_a.bedgraph";
        String pathB = TestUtils.DATA_DIR + "wig/toCombine_b.bedgraph";

        WiggleDataset dsA = (new WiggleParser(new ResourceLocator(pathA), genome)).parse();
        WiggleDataset dsB = new WiggleParser(new ResourceLocator(pathB), genome).parse();

        DataSource sourceA  = new DatasetDataSource("A", dsA, genome);
        DataSource sourceB  = new DatasetDataSource("B", dsB, genome);

        sourceA.setWindowFunction(WindowFunction.none);
        sourceB.setWindowFunction(WindowFunction.none);

        DataSourceTrack trackA = new DataSourceTrack(null, "A", "A", sourceA);
        DataSourceTrack trackB = new DataSourceTrack(null, "B", "B", sourceB);

        return new CombinedDataSource(trackA, trackB, operation);
    }

}
