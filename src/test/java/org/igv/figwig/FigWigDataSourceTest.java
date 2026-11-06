package org.igv.figwig;

import org.igv.AbstractHeadlessTest;
import org.igv.data.AbstractDataSource;
import org.igv.data.CompositeScore;
import org.igv.data.DataTile;
import org.igv.data.NamedScore;
import org.igv.feature.LocusScore;
import org.igv.feature.genome.Genome;
import org.igv.track.DataSourceTrack;
import org.igv.track.DataType;
import org.igv.track.Track;
import org.igv.track.TrackLoader;
import org.igv.track.WindowFunction;
import org.igv.util.ResourceLocator;
import org.json.JSONObject;
import org.junit.Test;

import java.lang.reflect.Field;
import java.util.*;

import static org.junit.Assert.*;

/**
 * figwig tracks through IGV: loading, attributes, releasing the file, and the summary of raw data,
 * which FigWigDataSource computes itself and must give exactly what AbstractDataSource gives.
 */
public class FigWigDataSourceTest extends AbstractHeadlessTest {

    static final String FILE = FigWigFileTest.DIR + "positions_zlib.fw";

    static FigWigDataSource source(Track track) throws Exception {
        Field f = DataSourceTrack.class.getDeclaredField("dataSource");
        f.setAccessible(true);
        return (FigWigDataSource) f.get(track);
    }

    @Test
    public void testLoadMakesATrackPerFigwigTrackWithItsAttributes() throws Exception {
        List<Track> tracks = new TrackLoader().load(new ResourceLocator(FILE), genome);
        FigWigFile fw = source(tracks.get(0)).getFile();
        assertEquals(fw.nTracks(), tracks.size());
        Map<String, Track> byName = new HashMap<>();
        for (Track t : tracks) {
            byName.put(t.getName(), t);
        }
        assertEquals("ATAC", byName.get("signal").getAttributeValue("assay"));
        assertEquals("K562", byName.get("signal").getAttributeValue("biosample"));
        assertEquals("17", byName.get("cell17").getAttributeValue("cell"));
        for (Track t : tracks) {
            t.unload();
        }
    }

    @Test
    public void testRemovingEveryTrackClosesTheFile() throws Exception {
        List<Track> tracks = new TrackLoader().load(new ResourceLocator(FILE), genome);
        FigWigFile fw = source(tracks.get(0)).getFile();
        for (int i = 0; i < tracks.size() - 1; i++) {
            tracks.get(i).unload();
        }
        assertFalse(fw.isClosed());
        tracks.get(tracks.size() - 1).unload();
        assertTrue(fw.isClosed());
    }

    /** IGV's own summary of the same raw and precomputed data: AbstractDataSource, not overridden. */
    static class Reference extends AbstractDataSource {
        final FigWigDataSource source;

        Reference(FigWigDataSource source, Genome genome) {
            super(genome);
            this.source = source;
        }

        @Override
        protected DataTile getRawData(String chr, int start, int end) {
            return source.getRawData(chr, start, end);
        }

        @Override
        protected List<LocusScore> getPrecomputedSummaryScores(String chr, int start, int end, int zoom) {
            source.setWindowFunction(getWindowFunction());
            return source.getPrecomputedSummaryScores(chr, start, end, zoom);
        }

        @Override
        public int getLongestFeature(String chr) {
            return 0;
        }

        public double getDataMax() {
            return source.getDataMax();
        }

        public double getDataMin() {
            return source.getDataMin();
        }

        public DataType getDataType() {
            return DataType.OTHER;
        }

        public void dispose() {
        }
    }

    static String describe(List<LocusScore> scores) throws Exception {
        Field data = CompositeScore.class.getDeclaredField("data");
        Field probes = CompositeScore.class.getDeclaredField("probes");
        Field probe = NamedScore.class.getDeclaredField("probe");
        data.setAccessible(true);
        probes.setAccessible(true);
        probe.setAccessible(true);
        StringBuilder sb = new StringBuilder();
        for (LocusScore s : scores) {
            sb.append(s.getClass().getSimpleName()).append(' ').append(s.getStart()).append(' ').append(s.getEnd())
                    .append(' ').append(Float.floatToRawIntBits(s.getScore()));
            if (s instanceof CompositeScore) {
                sb.append(" d=").append(Arrays.toString((float[]) data.get(s)))
                        .append(" p=").append(Arrays.toString((String[]) probes.get(s)));
            }
            if (s instanceof NamedScore) {
                sb.append(" n=").append(probe.get(s));
            }
            sb.append('\n');
        }
        return sb.toString();
    }

    @Test
    public void testSummaryIsExactlyAbstractDataSources() throws Exception {
        List<Track> tracks = new TrackLoader().load(new ResourceLocator(FILE), genome);
        int compared = 0;
        Set<String> some = new HashSet<>(Arrays.asList("cell0", "cell16", "cell32", "cell48", "site0", "site16"));
        for (Track t : tracks) {
            // The 15 named tracks hold every element type; of the cells and sites, a few.
            if ((t.getName().startsWith("cell") || t.getName().startsWith("site")) && !some.contains(t.getName())) {
                continue;
            }
            FigWigDataSource source = source(t);
            for (WindowFunction wf : new WindowFunction[]{WindowFunction.mean, WindowFunction.none}) {
                for (String chr : new String[]{"chr1", "chr2", "chr3"}) {
                    for (int zoom : new int[]{4, 10, 14, 18, 20, 22, 24}) {
                        for (int[] range : new int[][]{{0, 300}, {17, 251}, {60, 70}, {0, 1000000}}) {
                            source.setWindowFunction(wf);
                            List<LocusScore> fast = source.getSummaryScoresForRange(chr, range[0], range[1], zoom);
                            Reference ref = new Reference(source, genome);
                            ref.setWindowFunction(wf);
                            List<LocusScore> igv = ref.getSummaryScoresForRange(chr, range[0], range[1], zoom);
                            source.setWindowFunction(wf);
                            assertEquals(t.getName() + " " + wf + " " + chr + ":" + range[0] + "-" + range[1] + " zoom " + zoom,
                                    describe(igv), describe(fast));
                            compared += igv.size();
                        }
                    }
                }
            }
        }
        assertTrue("scores compared: " + compared, compared > 10000);
        for (Track t : tracks) {
            t.unload();
        }
    }

    @Test
    public void testChooserLoadsOnlyTheChosenTracks() throws Exception {
        try {
            FigWigTrackChooser.override = f -> new HashSet<>(Arrays.asList("signal", "ones", "cell17"));
            List<Track> tracks = new TrackLoader().load(new ResourceLocator(FILE), genome);
            Set<String> names = new HashSet<>();
            for (Track t : tracks) names.add(t.getName());
            assertEquals(new HashSet<>(Arrays.asList("signal", "ones", "cell17")), names);
            FigWigFile fw = source(tracks.get(0)).getFile();
            for (Track t : tracks) t.unload();
            assertTrue("removing the chosen tracks closes the file", fw.isClosed());

            FigWigTrackChooser.override = f -> new HashSet<>();   // cancelled
            assertEquals(0, new TrackLoader().load(new ResourceLocator(FILE), genome).size());
        } finally {
            FigWigTrackChooser.override = null;
        }
    }

    @Test
    public void testChooserShowsEveryTrackWithItsAttributes() throws Exception {
        FigWigFile fw = FigWigFile.open(FILE);
        try {
            List<String> headings = FigWigTrackChooser.headings(fw);
            assertEquals("name", headings.get(0));
            assertTrue(headings.containsAll(Arrays.asList("assay", "biosample", "cell")));
            List<org.igv.encode.FileRecord> records = FigWigTrackChooser.records(fw);
            assertEquals(fw.nTracks(), records.size());
            for (int t = 0; t < fw.nTracks(); t++) {
                assertEquals(fw.trackNames[t], records.get(t).getAttributeValue("name"));
            }
            assertTrue("small files load without asking: " + fw.nTracks(),
                    fw.nTracks() > FigWigTrackChooser.CHOOSER_MIN_TRACKS);
            assertNull("headless runs load every track", FigWigTrackChooser.choose(fw, "test"));
        } finally {
            fw.close();
        }
    }

    /**
     * In files where some groups keep a summary and others none: views of every track are exactly
     * AbstractDataSource's; the whole-genome view, the one view of these small chromosomes that uses a
     * summary, draws a track whose group keeps one and nothing for a track whose group keeps none.
     */
    @Test
    public void testGroupsWithoutASummaryUseRawData() throws Exception {
        for (String file : new String[]{"mixed.fw", "mixed_stored.fw", "mixed_merged.fw"}) {
            List<Track> tracks = new TrackLoader().load(new ResourceLocator(FigWigFileTest.DIR + file), genome);
            try {
                int with = 0, without = 0;
                for (Track t : tracks) {
                    if (!Arrays.asList("random", "signal", "cell0", "cell63", "site5", "bare").contains(t.getName())) {
                        continue;
                    }
                    FigWigDataSource source = source(t);
                    FigWigFile fw = source.getFile();
                    boolean has = fw.hasSummary(Arrays.asList(fw.trackNames).indexOf(t.getName()));
                    assertEquals(file + " " + t.getName(), !t.getName().startsWith("cell") && !t.getName().startsWith("site")
                            && !t.getName().equals("bare"), has);
                    source.setWindowFunction(WindowFunction.mean);
                    for (String chr : new String[]{"chr1", "chr2", "chr3"}) {
                        for (int zoom : new int[]{0, 1, 2, 4, 10}) {
                            List<LocusScore> fast = source.getSummaryScoresForRange(chr, 0, 300, zoom);
                            Reference ref = new Reference(source, genome);
                            ref.setWindowFunction(WindowFunction.mean);
                            assertEquals(file + " " + t.getName() + " " + chr + " zoom " + zoom,
                                    describe(ref.getSummaryScoresForRange(chr, 0, 300, zoom)), describe(fast));
                        }
                    }
                    List<LocusScore> whole = source.getSummaryScoresForRange(org.igv.Globals.CHR_ALL, 0,
                            (int) (genome.getWGLength() / 1000), 0);
                    if (has) {
                        assertFalse(file + " " + t.getName() + " whole genome", whole.isEmpty());
                        with++;
                    } else {
                        assertTrue(file + " " + t.getName() + " whole genome", whole.isEmpty());
                        without++;
                    }
                }
                assertTrue(file + ": " + with + " with a summary, " + without + " without", with > 0 && without > 0);
            } finally {
                for (Track t : tracks) t.unload();
            }
        }
    }

    @Test
    public void testWholeGenomeViewUsesASummaryLevelInGenomeCoordinates() throws Exception {
        List<Track> tracks = new TrackLoader().load(new ResourceLocator(FILE), genome);
        try {
            Track random = null;
            for (Track t : tracks) if (t.getName().equals("random")) random = t;
            FigWigDataSource source = source(random);
            source.setWindowFunction(WindowFunction.mean);
            assertEquals(org.igv.Globals.CHR_ALL, genome.getHomeChromosome());
            List<LocusScore> scores = source.getSummaryScoresForRange(org.igv.Globals.CHR_ALL, 0,
                    (int) (genome.getWGLength() / 1000), 0);
            assertFalse(scores.isEmpty());
            FigWigFile fw = source.getFile();
            int level = -1;
            for (int k = 0; k < fw.levels(); k++) {
                if (fw.levelBases(k) <= genome.getWGLength() / 1000.0) level = k;
            }
            int t = Arrays.asList(fw.trackNames).indexOf("random");
            FigWigFile.Runs chr1 = fw.levelMeans(t, level, fw.chromId("chr1"), 0, (int) fw.chromLengths[fw.chromId("chr1")]);
            assertEquals(genome.getGenomeCoordinate("chr1", chr1.starts[0]), scores.get(0).getStart());
            assertEquals((float) chr1.values[0], scores.get(0).getScore(), 0f);
            int chr3End = genome.getGenomeCoordinate("chr3", (int) fw.chromLengths[fw.chromId("chr3")]);
            for (LocusScore score : scores) {
                assertTrue(score.getStart() <= chr3End + 1);
            }
        } finally {
            for (Track t : tracks) t.unload();
        }
    }

    /**
     * A session saves each track with its file's path, format and the track's name; restoring loads the
     * file and takes the track of that name (JSONSessionReader), so a session keeps the tracks chosen.
     */
    @Test
    public void testSessionJsonFindsTheSameTrackAgain() throws Exception {
        List<Track> tracks = new TrackLoader().load(new ResourceLocator(FILE), genome);
        Track signal = null;
        for (Track t : tracks) if (t.getName().equals("signal")) signal = t;
        ((DataSourceTrack) signal).setWindowFunction(WindowFunction.none);
        JSONObject json = new JSONObject();
        signal.marshalJSON(json);
        for (Track t : tracks) t.unload();

        assertEquals(new java.io.File(FILE).getPath(), new java.io.File(json.getString("path")).getPath());
        assertEquals("fw", json.getString("format"));
        assertEquals("signal", json.getString("name"));

        ResourceLocator locator = new ResourceLocator(json.getString("path"));
        locator.setFormat(json.getString("format"));
        List<Track> again = new TrackLoader().load(locator, genome);
        Track match = null;
        for (Track t : again) if (t.getName().equals(json.getString("name"))) match = t;
        assertNotNull(match);
        match.unmarshalJSON(json);
        assertEquals(WindowFunction.none, ((DataSourceTrack) match).getWindowFunction());
        assertEquals(json.getString("id"), match.getId());
        for (Track t : again) t.unload();
    }
}
