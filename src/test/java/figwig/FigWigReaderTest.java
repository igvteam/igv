package figwig;

import org.junit.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.*;

/** FigWigReader against the values the figwig package reads from the same files (see FigWigFileTest). */
public class FigWigReaderTest {

    /** Windows in any order, overlapping and running past a chromosome's end, from a reader of a few tracks. */
    @Test
    public void testWindowsHoldWhatTheFigwigPackageReads() throws Exception {
        for (String file : new String[]{"positions_zlib.fw", "tracks_zlib.fw", "extremes.fw"}) {
            FigWigFileTest.Expected e = FigWigFileTest.expected(file);
            FigWigFile all = FigWigFile.open(FigWigFileTest.DIR + file);
            String[] names = {all.trackName(all.nTracks() - 1), all.trackName(0), all.trackName(1)};
            int[] fileIndex = {all.nTracks() - 1, 0, 1};
            try (FigWigReader fw = new FigWigReader(all, names)) {
                String last = all.chromName(all.nChroms() - 1);
                int lastLength = (int) all.chromLength(all.nChroms() - 1);
                String[] chroms = {last, all.chromName(0), all.chromName(0), last};
                int[] starts = {Math.max(0, lastLength - 40), 70, 0, 5};
                int width = 100;
                double[][][] y = fw.read(chroms, starts, width);
                assertEquals(4, y.length);
                assertEquals(3, y[0].length);
                for (int j = 0; j < chroms.length; j++) {
                    long length = all.chromLength(all.chromId(chroms[j]));
                    for (int i = 0; i < names.length; i++) {
                        Map<Integer, Long> want = e.of(fileIndex[i], chroms[j]);
                        for (int p = 0; p < width; p++) {
                            int at = starts[j] + p;
                            String where = file + " " + names[i] + " " + chroms[j] + ":" + at;
                            if (at >= length) {
                                assertTrue(where, Double.isNaN(y[j][i][p]));
                            } else {
                                assertEquals(where, (long) want.getOrDefault(at, Double.doubleToRawLongBits(0.0)),
                                        FigWigFileTest.bits(y[j][i][p]));
                            }
                        }
                    }
                }
                // one track, and tracks in another order
                assertArrayEquals(y[1][2], fw.read(chroms, starts, width, 2)[1], 0.0);
                assertArrayEquals(y[3][0], fw.read(chroms, starts, width, new int[]{1, 0})[3][1], 0.0);
            }
            assertTrue("closing the reader closes the file", all.isClosed());
        }
    }

    /** The sums the file keeps are the sums of its values: per tile, and per bin of each zoom level. */
    @Test
    public void testTileSumsAreTheSumsOfTheValues() throws Exception {
        try (FigWigReader fw = new FigWigReader(FigWigFileTest.DIR + "positions_zlib.fw")) {
            int tile = fw.tile();
            int[] zooms = fw.zooms();
            assertTrue(zooms.length > 0);
            int checked = 0;
            for (int t : fw.select("dtype", "uint8")) {
                for (Map.Entry<String, Long> chrom : fw.chromSizes().entrySet()) {
                    int length = (int) (long) chrom.getValue();
                    double[] values = fw.read(new String[]{chrom.getKey()}, new int[]{0}, length, t)[0];
                    double[] sums = fw.tileSums(t, chrom.getKey(), 0);
                    assertEquals((length + tile - 1) / tile, sums.length);
                    for (int k = 0; k < sums.length; k++) {
                        double want = 0;
                        for (int p = k * tile; p < Math.min(length, (k + 1) * tile); p++) {
                            want += values[p];
                        }
                        assertEquals(chrom.getKey() + " tile " + k, want, sums[k], 0.0);
                        checked += want != 0 ? 1 : 0;
                    }
                    for (int z = 0; z < zooms.length; z++) {
                        double[] zoomed = fw.tileSums(t, chrom.getKey(), z + 1);
                        assertEquals((sums.length + zooms[z] - 1) / zooms[z], zoomed.length);
                        for (int b = 0; b < zoomed.length; b++) {
                            double want = 0;
                            for (int k = b * zooms[z]; k < Math.min(sums.length, (b + 1) * zooms[z]); k++) {
                                want += sums[k];
                            }
                            assertEquals(chrom.getKey() + " level " + (z + 1) + " bin " + b, want, zoomed[b], 0.0);
                        }
                    }
                }
            }
            assertTrue(checked > 10);
        }
        try (FigWigReader cells = new FigWigReader(FigWigFileTest.DIR + "mixed.fw", "cell0")) {
            cells.tileSums(0, "chr1", 0);
            fail("a group that keeps no summary gave sums");
        } catch (IllegalStateException expected) {
            // no sums are kept for this track
        }
    }

    @Test
    public void testSelectMatchesTracksByTheirAttributes() throws Exception {
        try (FigWigReader fw = new FigWigReader(FigWigFileTest.DIR + "positions_zlib.fw")) {
            List<Map<String, String>> tracks = fw.tracks();
            int[] atac = fw.select("assay", "ATAC", "biosample", "K562");
            assertTrue(atac.length > 0 && atac.length < tracks.size());
            boolean signal = false;
            for (int t : atac) {
                assertEquals("ATAC", tracks.get(t).get("assay"));
                signal |= tracks.get(t).get("name").equals("signal");
            }
            assertTrue(signal);
            assertEquals(0, fw.select("assay", "no such assay").length);
            assertEquals("[" + Arrays.stream(fw.select("cell", "17")).mapToObj(t -> tracks.get(t).get("name"))
                    .reduce("", String::concat) + "]", "[cell17]");
        }
    }
}
