package org.igv.figwig;

import org.igv.AbstractHeadlessTest;
import org.igv.feature.LocusScore;
import org.igv.track.Track;
import org.igv.track.TrackLoader;
import org.igv.track.WindowFunction;
import org.igv.util.ResourceLocator;
import org.junit.Test;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.*;

import static org.junit.Assert.*;

/**
 * Cases a user meets outside the common path: every window function, chromosome names that differ from
 * the genome's, requests reaching past a chromosome, extreme values, and reads racing a file's closing.
 */
public class FigWigEdgeCaseTest extends AbstractHeadlessTest {

    /** Every window function IGV offers anywhere, at raw zooms, is exactly AbstractDataSource's. */
    @Test
    public void testEveryWindowFunctionIsAbstractDataSources() throws Exception {
        for (String file : new String[]{"positions_zlib.fw", "extremes.fw"}) {
            List<Track> tracks = new TrackLoader().load(new ResourceLocator(FigWigFileTest.DIR + file), genome);
            try {
                int compared = 0;
                for (Track t : tracks) {
                    if (t.getName().startsWith("cell") || t.getName().startsWith("site")) {
                        continue;
                    }
                    FigWigDataSource source = FigWigDataSourceTest.source(t);
                    for (WindowFunction wf : WindowFunction.values()) {
                        for (String chr : new String[]{"chr1", "chr2"}) {
                            for (int zoom : new int[]{0, 6, 12, 20}) {
                                source.setWindowFunction(wf);
                                List<LocusScore> fast = source.getSummaryScoresForRange(chr, 0, 400, zoom);
                                FigWigDataSourceTest.Reference ref = new FigWigDataSourceTest.Reference(source, genome);
                                ref.setWindowFunction(wf);
                                List<LocusScore> igv = ref.getSummaryScoresForRange(chr, 0, 400, zoom);
                                source.setWindowFunction(wf);
                                assertEquals(file + " " + t.getName() + " " + wf + " " + chr + " zoom " + zoom,
                                        FigWigDataSourceTest.describe(igv), FigWigDataSourceTest.describe(fast));
                                compared += igv.size();
                            }
                        }
                    }
                }
                assertTrue(compared > 1000);
            } finally {
                for (Track t : tracks) t.unload();
            }
        }
    }

    /** A file named in Ensembl style (1, 2, MT) draws under the genome's names for them. */
    @Test
    public void testChromosomeNamesFollowTheGenomesAliases() throws Exception {
        List<Track> tracks = new TrackLoader().load(new ResourceLocator(FigWigFileTest.DIR + "ensembl.fw"), genome);
        try {
            FigWigDataSource source = FigWigDataSourceTest.source(tracks.get(1));
            FigWigFile fw = source.getFile();
            String chr1 = genome.getCanonicalChrName("1");
            assertNotEquals("1", chr1);
            source.setWindowFunction(WindowFunction.none);
            List<LocusScore> scores = source.getSummaryScoresForRange(chr1, 0, 300, 20);
            FigWigFile.Runs runs = fw.read(1, fw.chromId("1"), 0, 300);
            assertFalse(scores.isEmpty());
            assertEquals(runs.n, scores.size());
            for (int i = 0; i < runs.n; i++) {
                assertEquals(runs.starts[i], scores.get(i).getStart());
                assertEquals((float) runs.values[i], scores.get(i).getScore(), 0f);
            }
            assertTrue(source.getSummaryScoresForRange("chrNotInTheFile", 0, 300, 20).isEmpty());
        } finally {
            for (Track t : tracks) t.unload();
        }
    }

    /** Requests starting before 0, ending past the chromosome, empty, or wholly past it. */
    @Test
    public void testRangesReachingPastAChromosome() throws Exception {
        for (String file : new String[]{"positions_zlib.fw", "extremes.fw", "tracks_zlib.fw"}) {
            FigWigFile fw = FigWigFile.open(FigWigFileTest.DIR + file);
            FigWigFileTest.Expected e = FigWigFileTest.expected(file);
            try {
                for (int t = 0; t < fw.nTracks(); t++) {
                    for (int c = 0; c < fw.chromNames.length; c++) {
                        int len = (int) fw.chromLengths[c];
                        for (int[] r : new int[][]{{-50, 10}, {len - 5, len + 100}, {-10, len + 10}, {len, len + 10},
                                {len + 5, len + 50}, {7, 7}, {9, 3}}) {
                            FigWigFile.Runs runs = fw.read(t, c, r[0], r[1]);
                            int from = Math.max(0, r[0]), to = Math.min(len, r[1]);
                            Map<Integer, Long> want = e.of(t, fw.chromNames[c]);
                            int expectedRuns = 0;
                            for (int i = 0; i < runs.n; i++) {
                                assertTrue(runs.starts[i] >= from && runs.ends[i] <= to && runs.starts[i] < runs.ends[i]);
                                for (int p = runs.starts[i]; p < runs.ends[i]; p++) {
                                    assertEquals(file + " " + t + " " + fw.chromNames[c] + ":" + p,
                                            (long) want.get(p), FigWigFileTest.bits(runs.values[i]));
                                }
                            }
                            for (int p = from; p < to; p++) {
                                if (want.containsKey(p)) expectedRuns++;
                            }
                            int covered = 0;
                            for (int i = 0; i < runs.n; i++) covered += runs.ends[i] - runs.starts[i];
                            assertEquals(file + " " + t + " " + fw.chromNames[c] + " " + Arrays.toString(r), expectedRuns, covered);
                        }
                    }
                }
            } finally {
                fw.close();
            }
        }
    }

    /** Reads that race the file's closing either return the right values or fail with IOException. */
    @Test(timeout = 120000)
    public void testReadsRacingCloseFailCleanly() throws Exception {
        FigWigFileTest.Expected e = FigWigFileTest.expected("cells.fw");
        for (int round = 0; round < 20; round++) {
            FigWigFile fw = FigWigFile.open(FigWigFileTest.DIR + "cells.fw");
            for (int t = 0; t < fw.nTracks(); t++) fw.retain(t);
            ExecutorService pool = Executors.newFixedThreadPool(8);
            List<Future<String>> results = new ArrayList<>();
            Random rng = new Random(round);
            for (int i = 0; i < 400; i++) {
                final int t = rng.nextInt(fw.nTracks()), c = rng.nextInt(fw.chromNames.length);
                final int len = (int) fw.chromLengths[c];
                results.add(pool.submit(() -> {
                    try {
                        double[] v = fw.readDense(t, c, 0, len);
                        Map<Integer, Long> want = e.of(t, fw.chromNames[c]);
                        for (int p = 0; p < len; p++) {
                            if (FigWigFileTest.bits(v[p]) != want.getOrDefault(p, Double.doubleToRawLongBits(0.0))) {
                                return "wrong value at track " + t + " " + fw.chromNames[c] + ":" + p;
                            }
                        }
                        return null;
                    } catch (IOException | java.io.UncheckedIOException ok) {
                        return null;
                    }
                }));
            }
            Thread.sleep(rng.nextInt(5));
            for (int t = 0; t < fw.nTracks(); t++) fw.release(t);
            assertTrue(fw.isClosed());
            for (Future<String> f : results) {
                try {
                    assertNull(f.get(), f.get());
                } catch (ExecutionException x) {
                    fail("a read racing close threw " + x.getCause());
                }
            }
            pool.shutdown();
        }
    }
}
