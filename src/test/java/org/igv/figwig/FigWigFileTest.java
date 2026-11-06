package org.igv.figwig;

import org.igv.util.TestUtils;
import org.junit.Test;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.*;
import java.util.zip.GZIPInputStream;

import static org.junit.Assert.*;

/**
 * The reader against the values the figwig package reads from the same files.  The files in
 * test/data/figwig were written by the package, and each FILE.expected.tsv.gz lists every nonzero base
 * of every track and every nonzero bin of every summary and zoom level as the package reads them.
 * Between them the files hold every chunk kind and element type, both orders and both codecs, and in
 * mixed*.fw groups that keep no summary (listed as nosummary, with no levels).
 */
public class FigWigFileTest {

    static final String DIR = TestUtils.DATA_DIR + "figwig/";
    static final String[] FILES = {"positions_zlib.fw", "positions_stored.fw", "tracks_zlib.fw", "tracks_stored.fw",
            "entries.fw", "long.fw", "cells.fw", "extremes.fw", "ensembl.fw", "indexed.fw", "mixed.fw",
            "mixed_stored.fw", "mixed_merged.fw"};

    /** Expected values: per track and chromosome, position -> float64 bits (NaN as Long.MIN_VALUE). */
    static class Expected {
        final Map<String, Map<Integer, Long>> values = new HashMap<>();
        final Map<String, Long> levels = new HashMap<>();   // track:level:chrom:binStart -> bits
        final Set<Integer> noSummary = new HashSet<>();
        int maxLevel = -1;

        Map<Integer, Long> of(int track, String chrom) {
            return values.getOrDefault(track + ":" + chrom, Collections.emptyMap());
        }
    }

    static final long NAN = Long.MIN_VALUE;

    static Expected expected(String file) throws IOException {
        Expected e = new Expected();
        String path = DIR + file.replace(".fw", ".expected.tsv.gz");
        try (BufferedReader r = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(path))))) {
            String line;
            while ((line = r.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] f = line.split("\t");
                long bits = f[f.length - 1].equals("nan") ? NAN : Long.parseUnsignedLong(f[f.length - 1], 16);
                if (f[0].equals("nosummary")) {
                    e.noSummary.add(Integer.parseInt(f[1]));
                    continue;
                }
                if (f[0].equals("value")) {
                    e.values.computeIfAbsent(f[1] + ":" + f[2], k -> new HashMap<>()).put(Integer.parseInt(f[3]), bits);
                } else {
                    e.levels.put(f[1] + ":" + f[2] + ":" + f[3] + ":" + f[4], bits);
                    e.maxLevel = Math.max(e.maxLevel, Integer.parseInt(f[2]));
                }
            }
        }
        return e;
    }

    static long bits(double v) {
        return Double.isNaN(v) ? NAN : Double.doubleToRawLongBits(v);
    }

    static void assertWindow(Expected e, FigWigFile fw, int t, int c, int start, int end, String file) throws IOException {
        double[] got = fw.readDense(t, c, start, end);
        Map<Integer, Long> want = e.of(t, fw.chromNames[c]);
        for (int p = start; p < end; p++) {
            long w = want.getOrDefault(p, Double.doubleToRawLongBits(0.0));
            assertEquals(file + " track " + t + " " + fw.chromNames[c] + ":" + p, w, bits(got[p - start]));
        }
    }

    @Test
    public void testEveryBaseMatchesTheFigwigPackage() throws Exception {
        for (String file : FILES) {
            Expected e = expected(file);
            FigWigFile fw = FigWigFile.open(DIR + file);
            try {
                for (int t = 0; t < fw.nTracks(); t++) {
                    for (int c = 0; c < fw.chromNames.length; c++) {
                        assertWindow(e, fw, t, c, 0, (int) fw.chromLengths[c], file);
                    }
                }
            } finally {
                fw.close();
            }
        }
    }

    @Test
    public void testWindowsMatchTheFigwigPackage() throws Exception {
        Random rng = new Random(1);
        for (String file : FILES) {
            Expected e = expected(file);
            FigWigFile fw = FigWigFile.open(DIR + file);
            try {
                for (int i = 0; i < 400; i++) {
                    int t = rng.nextInt(fw.nTracks());
                    int c = rng.nextInt(fw.chromNames.length);
                    int len = (int) fw.chromLengths[c];
                    int start = rng.nextInt(len);
                    int end = start + 1 + rng.nextInt(len - start);
                    assertWindow(e, fw, t, c, start, end, file);
                }
            } finally {
                fw.close();
            }
        }
    }

    /**
     * Reads of many tracks of one group, as IGV makes them for a view of many cells: with every track in
     * use, the first read of a range decodes them all in one pass, and every track's values are still
     * the package's, for ranges long enough to batch and short ones read alone.
     */
    @Test
    public void testManyTracksOfAGroupReadInOnePass() throws Exception {
        String file = "cells.fw";
        Expected e = expected(file);
        FigWigFile fw = FigWigFile.open(DIR + file);
        ExecutorService pool = Executors.newFixedThreadPool(5);
        try {
            for (int t = 0; t < fw.nTracks(); t++) {
                fw.retain(t);
            }
            Random rng = new Random(3);
            int batched = 0;
            for (int i = 0; i < 60; i++) {
                int c = i % 2 == 0 ? 0 : rng.nextInt(fw.chromNames.length);
                int len = (int) fw.chromLengths[c];
                int start = i < 4 ? 0 : rng.nextInt(len);
                int end = i < 4 ? len : start + 1 + rng.nextInt(len - start);
                fw.clearCache();
                List<Future<double[]>> got = new ArrayList<>();
                for (int t = 0; t < fw.nTracks(); t++) {
                    final int track = t;
                    got.add(pool.submit(() -> fw.readDense(track, c, start, end)));
                }
                for (int t = 0; t < fw.nTracks(); t++) {
                    double[] v = got.get(t).get();
                    Map<Integer, Long> want = e.of(t, fw.chromNames[c]);
                    for (int p = start; p < end; p++) {
                        assertEquals(file + " track " + t + " " + fw.chromNames[c] + ":" + p + " in " + start + "-" + end,
                                (long) want.getOrDefault(p, Double.doubleToRawLongBits(0.0)), bits(v[p - start]));
                    }
                }
                batched += fw.batchesKept();
            }
            assertTrue("reads decoded in one pass", batched > 0);
        } finally {
            pool.shutdownNow();
            fw.close();
        }
    }

    /**
     * A file with a regional summary index: a chromosome's levels come from its windows alone and match
     * the package's; an index whose copy of the header is not the file's (stale) is ignored.
     */
    @Test
    public void testRegionalSummaryIndex() throws Exception {
        String file = "indexed.fw";
        Expected e = expected(file);
        FigWigFile fw = FigWigFile.open(DIR + file);
        try {
            for (int k = 0; k < fw.levels(); k++) {
                FigWigFile.Runs runs = fw.levelMeans(3, k, 2, 0, (int) fw.chromLengths[2]);
                for (int i = 0; i < runs.n; i++) {
                    for (long b = runs.starts[i]; b < runs.ends[i]; b += fw.levelBases(k)) {
                        assertEquals((long) e.levels.get("3:" + k + ":" + fw.chromNames[2] + ":" + b), bits(runs.values[i]));
                    }
                }
            }
            assertTrue("windows read", fw.summaryWindowsKept() > 0);
        } finally {
            fw.close();
        }
        // The stale copy: one byte of the root's header snapshot changed, the root's CRC-32 made again.
        byte[] b = Files.readAllBytes(new File(DIR + file).toPath());
        String meta = new String(b, java.nio.charset.StandardCharsets.ISO_8859_1);
        int key = meta.indexOf("\"offset\": ", meta.indexOf("org.figwig.summary-restarts")) + "\"offset\": ".length();
        int root = Integer.parseInt(meta.substring(key, meta.indexOf(',', key)).trim());
        b[root + 64 + 100] ^= 4;
        Arrays.fill(b, root + 12, root + 16, (byte) 0);
        java.util.zip.CRC32 crc = new java.util.zip.CRC32();
        crc.update(b, root, 256);
        java.nio.ByteBuffer.wrap(b, root + 12, 4).order(java.nio.ByteOrder.LITTLE_ENDIAN).putInt((int) crc.getValue());
        Path stale = Files.createTempFile("figwig-stale", ".fw");
        try {
            Files.write(stale, b);
            FigWigFile s = FigWigFile.open(stale.toString());
            try {
                for (int k = 0; k < s.levels(); k++) {
                    FigWigFile.Runs runs = s.levelMeans(3, k, 2, 0, (int) s.chromLengths[2]);
                    for (int i = 0; i < runs.n; i++) {
                        for (long x = runs.starts[i]; x < runs.ends[i]; x += s.levelBases(k)) {
                            assertEquals((long) e.levels.get("3:" + k + ":" + s.chromNames[2] + ":" + x), bits(runs.values[i]));
                        }
                    }
                }
                assertEquals("a stale index is not used", 0, s.summaryWindowsKept());
            } finally {
                s.close();
            }
        } finally {
            Files.delete(stale);
        }
    }

    @Test
    public void testSummaryAndZoomLevelsMatchTheFigwigPackage() throws Exception {
        for (String file : FILES) {
            Expected e = expected(file);
            FigWigFile fw = FigWigFile.open(DIR + file);
            try {
                int checked = 0;
                assertTrue(file + ": levels the package reads", e.maxLevel < fw.levels());
                for (int t = 0; t < fw.nTracks(); t++) {
                    assertEquals(file + " track " + t + " keeps a summary", !e.noSummary.contains(t), fw.hasSummary(t));
                    if (e.noSummary.contains(t)) {
                        for (int k = 0; k < fw.levels(); k++) {
                            assertNull(fw.levelMeans(t, k, 0, 0, (int) fw.chromLengths[0]));
                        }
                        checked++;
                        continue;
                    }
                    for (int k = 0; k < fw.levels(); k++) {
                        long bases = fw.levelBases(k);
                        for (int c = 0; c < fw.chromNames.length; c++) {
                            FigWigFile.Runs runs = fw.levelMeans(t, k, c, 0, (int) fw.chromLengths[c]);
                            Map<Long, Long> got = new HashMap<>();
                            for (int i = 0; i < runs.n; i++) {
                                // Runs merge neighbouring bins of one mean.
                                for (long b = runs.starts[i]; b < runs.ends[i]; b += bases) {
                                    got.put(b, bits(runs.values[i]));
                                }
                            }
                            for (long b = 0; b < fw.chromLengths[c]; b += bases) {
                                Long want = e.levels.get(t + ":" + k + ":" + fw.chromNames[c] + ":" + b);
                                assertEquals(file + " track " + t + " level " + k + " " + fw.chromNames[c] + ":" + b,
                                        want == null ? Double.doubleToRawLongBits(0.0) : want,
                                        (long) got.getOrDefault(b, Double.doubleToRawLongBits(0.0)));
                                checked++;
                            }
                        }
                    }
                }
                assertTrue(checked > 0);
            } finally {
                fw.close();
            }
        }
    }

    /** IGV loads the tracks of a view on several threads, which share the chunks of a group. */
    @Test(timeout = 120000)
    public void testConcurrentReadsMatchTheFigwigPackage() throws Exception {
        for (String file : FILES) {
            Expected e = expected(file);
            FigWigFile fw = FigWigFile.open(DIR + file);
            ExecutorService pool = Executors.newFixedThreadPool(8);
            try {
                for (int round = 0; round < 20; round++) {
                    fw.clearCache();
                    List<Future<?>> futures = new ArrayList<>();
                    for (int t = 0; t < fw.nTracks(); t++) {
                        final int track = t;
                        futures.add(pool.submit(() -> {
                            for (int c = 0; c < fw.chromNames.length; c++) {
                                assertWindow(e, fw, track, c, 0, (int) fw.chromLengths[c], file);
                            }
                            return null;
                        }));
                    }
                    for (Future<?> f : futures) {
                        f.get();
                    }
                }
            } finally {
                pool.shutdown();
                fw.close();
            }
        }
    }

    /**
     * Corrupted copies are refused with an IOException or IllegalStateException, or read; nothing else
     * escapes and nothing hangs.
     */
    @Test(timeout = 300000)
    public void testCorruptFilesFailCleanly() throws Exception {
        Path tmp = Files.createTempDirectory("figwig-corrupt");
        try {
            for (String file : new String[]{"positions_zlib.fw", "positions_stored.fw", "tracks_zlib.fw", "tracks_stored.fw",
                    "indexed.fw"}) {
                byte[] original = Files.readAllBytes(new File(DIR + file).toPath());
                List<byte[]> copies = new ArrayList<>();
                for (int at = 0; at < original.length; at += Math.max(1, original.length / 150)) {
                    byte[] b = original.clone();
                    b[at] ^= (byte) 0x5a;
                    copies.add(b);
                }
                for (int len : new int[]{0, 50, 127, 128, original.length / 3, original.length / 2, original.length - 1}) {
                    copies.add(Arrays.copyOf(original, len));
                }
                byte[] header = original.clone();
                header[104] ^= 1;   // the header's CRC-32
                copies.add(header);
                int refused = 0;
                List<String> unexpected = new ArrayList<>();
                for (int i = 0; i < copies.size(); i++) {
                    Path p = tmp.resolve(i + ".fw");
                    Files.write(p, copies.get(i));
                    try {
                        FigWigFile fw = FigWigFile.open(p.toString());
                        try {
                            for (int t = 0; t < fw.nTracks(); t++) {
                                for (int c = 0; c < fw.chromNames.length; c++) {
                                    fw.readDense(t, c, 0, (int) fw.chromLengths[c]);
                                }
                                for (int k = 0; k < fw.levels(); k++) {
                                    fw.levelMeans(t, k, 0, 0, (int) fw.chromLengths[0]);
                                }
                            }
                        } finally {
                            fw.close();
                        }
                    } catch (IOException | IllegalStateException | UncheckedIOException ex) {
                        refused++;
                    } catch (RuntimeException ex) {
                        unexpected.add("copy " + i + ": " + ex);
                    } finally {
                        Files.delete(p);
                    }
                }
                assertEquals(file + ": exceptions other than IOException and IllegalStateException",
                        Collections.emptyList(), unexpected);
                assertTrue(file + ": the truncated copies and the header CRC are refused", refused >= 8);
            }
        } finally {
            Files.deleteIfExists(tmp);
        }
    }

    @Test
    public void testARunAcrossManyTilesIsOneRun() throws Exception {
        // A read of all 200 tiles is split into parts loaded in parallel; their runs must join into one.
        FigWigFile fw = FigWigFile.open(DIR + "long.fw");
        try {
            int c = fw.chromId("chr1");
            FigWigFile.Runs runs = fw.read(0, c, 0, (int) fw.chromLengths[c]);
            assertEquals(1, runs.n);
            assertEquals(100, runs.starts[0]);
            assertEquals(100 + 140 * 64, runs.ends[0]);
            assertEquals(2.5, runs.values[0], 0);
        } finally {
            fw.close();
        }
    }

    @Test
    public void testTheLastReleaseClosesTheFile() throws Exception {
        FigWigFile fw = FigWigFile.open(DIR + "positions_zlib.fw");
        fw.retain();
        fw.retain();
        fw.release();
        assertFalse(fw.isClosed());
        fw.readDense(0, 0, 0, 10);
        fw.release();
        assertTrue(fw.isClosed());
        try {
            fw.clearCache();
            fw.readDense(1, 0, 0, 10);
            fail("a closed file was read");
        } catch (IOException expected) {
            // closed
        }
    }
}
