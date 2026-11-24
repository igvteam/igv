package org.igv.figwig;

import figwig.FigWigFile;
import org.igv.AbstractHeadlessTest;
import org.igv.feature.LocusScore;
import org.igv.track.Track;
import org.igv.track.TrackLoader;
import org.igv.track.WindowFunction;
import org.igv.util.ResourceLocator;
import org.junit.Test;

import java.util.*;

import static org.junit.Assert.*;

/**
 * Cases a user meets outside the common path: every window function, extreme values, and chromosome
 * names that differ from the genome's.
 */
public class FigWigEdgeCaseTest extends AbstractHeadlessTest {

    /** Every window function IGV offers anywhere, at raw zooms, is exactly AbstractDataSource's. */
    @Test
    public void testEveryWindowFunctionIsAbstractDataSources() throws Exception {
        for (String file : new String[]{"positions_zlib.fw", "extremes.fw"}) {
            List<Track> tracks = new TrackLoader().load(new ResourceLocator(FigWigDataSourceTest.DIR + file), genome);
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
        List<Track> tracks = new TrackLoader().load(new ResourceLocator(FigWigDataSourceTest.DIR + "ensembl.fw"), genome);
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
}
