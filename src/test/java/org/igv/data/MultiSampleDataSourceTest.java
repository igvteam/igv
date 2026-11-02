package org.igv.data;

import org.igv.AbstractHeadlessTest;
import org.igv.Globals;
import org.igv.feature.LocusScore;
import org.igv.tdf.TDFDataSource;
import org.igv.tdf.TDFMultiSampleSource;
import org.igv.tdf.TDFReader;
import org.igv.track.WindowFunction;
import org.igv.util.ResourceLocator;
import org.igv.util.TestUtils;
import org.junit.Test;

import java.util.List;
import java.util.Map;

import static org.junit.Assert.*;

/**
 * A multi-sample source fetches a region once for all samples.  Each sample's scores must equal what the
 * single-sample sources produce for that sample alone.
 */
public class MultiSampleDataSourceTest extends AbstractHeadlessTest {

    @Test
    public void testDatasetMatchesPerSampleSource() throws Exception {

        IGVDataset ds = new IGVDataset(new ResourceLocator(TestUtils.DATA_DIR + "igv/MIP_44.cn"), genome);
        String[] samples = ds.getSampleNames();
        assertEquals(44, samples.length);
        DatasetMultiSampleSource source = new DatasetMultiSampleSource(ds, genome);
        assertEquals(samples.length, source.getSampleNames().size());

        int chr1Length = genome.getChromosome("chr1").getLength();
        int[][] regions = {{0, chr1Length / 8, 4}, {chr1Length / 3, chr1Length / 3 + 2000000, 10}};
        for (WindowFunction wf : new WindowFunction[]{WindowFunction.mean, WindowFunction.max, WindowFunction.none}) {
            source.setWindowFunction(wf);
            for (int[] region : regions) {
                Map<String, List<LocusScore>> all = source.getFeatures("chr1", region[0], region[1], region[2]);
                assertEquals(samples.length, all.size());
                for (String sample : new String[]{samples[0], samples[7], samples[43]}) {
                    DatasetDataSource single = new DatasetDataSource(sample, ds, genome);
                    single.setWindowFunction(wf);
                    List<LocusScore> expected = single.getSummaryScoresForRange("chr1", region[0], region[1], region[2]);
                    assertFalse(expected.isEmpty());
                    assertScoresEqual(wf + " " + sample, expected, all.get(sample));
                }
            }
        }

        // Whole genome: one score per genome summary point, at the summary's own resolution
        source.setWindowFunction(WindowFunction.mean);
        Map<String, List<LocusScore>> wg = source.getFeatures(Globals.CHR_ALL, 0, Integer.MAX_VALUE, 0);
        GenomeSummaryData summary = ds.getGenomeSummary();
        int[] starts = summary.getLocations();
        int[] ends = summary.getEndLocations();
        float[] values = summary.getData(samples[3]);
        List<LocusScore> actual = wg.get(samples[3]);
        assertEquals(starts.length, actual.size());
        for (int i = 0; i < starts.length; i++) {
            assertEquals(starts[i], actual.get(i).getStart());
            assertEquals(ends[i], actual.get(i).getEnd());
            assertEquals(values[i], actual.get(i).getScore(), 0f);
            assertTrue(ends[i] > starts[i]);
        }
    }

    @Test
    public void testTDFMatchesPerTrackSource() throws Exception {

        ResourceLocator locator = new ResourceLocator(TestUtils.DATA_DIR + "igv/recombRate.igv.txt.tdf");
        TDFReader reader = TDFReader.getReader(locator);
        String[] tracks = reader.getSampleNames();
        assertEquals(9, tracks.length);
        TDFMultiSampleSource source = new TDFMultiSampleSource(reader, genome);

        // Precomputed zoom level
        int chr1Length = genome.getChromosome("chr1").getLength();
        int start = chr1Length / 4;
        int end = chr1Length / 2;
        Map<String, List<LocusScore>> all = source.getFeatures("chr1", start, end, 3);
        for (int t = 0; t < tracks.length; t++) {
            TDFDataSource single = new TDFDataSource(TDFReader.getReader(locator), t, tracks[t], genome);
            List<LocusScore> expected = single.getSummaryScoresForRange("chr1", start, end, 3);
            assertFalse(expected.isEmpty());
            assertScoresEqual(tracks[t], expected, all.get(tracks[t]));
        }

        // Whole genome
        Map<String, List<LocusScore>> wg = source.getFeatures(Globals.CHR_ALL, 0, Integer.MAX_VALUE, 0);
        TDFDataSource single = new TDFDataSource(TDFReader.getReader(locator), 2, tracks[2], genome);
        List<LocusScore> expected = single.getSummaryScoresForRange(Globals.CHR_ALL, 0, Integer.MAX_VALUE, 0);
        assertFalse(expected.isEmpty());
        assertScoresEqual("whole genome", expected, wg.get(tracks[2]));

        // Raw data (zoom beyond the precomputed levels), no window function: same features and values per track
        source.setWindowFunction(WindowFunction.none);
        single.setWindowFunction(WindowFunction.none);
        all = source.getFeatures("chr1", start, start + 5000000, 9);
        expected = single.getSummaryScoresForRange("chr1", start, start + 5000000, 9);
        List<LocusScore> actual = all.get(tracks[2]);
        assertFalse(expected.isEmpty());
        assertEquals(expected.size(), actual.size());
        for (int i = 0; i < expected.size(); i++) {
            assertEquals(expected.get(i).getStart(), actual.get(i).getStart());
            assertEquals(expected.get(i).getScore(), actual.get(i).getScore(), 0);
        }
    }

    private static void assertScoresEqual(String message, List<LocusScore> expected, List<LocusScore> actual) {
        assertNotNull(message, actual);
        assertEquals(message + " count", expected.size(), actual.size());
        for (int i = 0; i < expected.size(); i++) {
            LocusScore e = expected.get(i);
            LocusScore a = actual.get(i);
            assertEquals(message + " start " + i, e.getStart(), a.getStart());
            assertEquals(message + " end " + i, e.getEnd(), a.getEnd());
            assertEquals(message + " score " + i, e.getScore(), a.getScore(), 0);
        }
    }
}
