package org.igv.figwig;

import org.igv.Globals;
import org.igv.data.AbstractDataSource;
import org.igv.data.BasicScore;
import org.igv.data.CompositeScore;
import org.igv.data.DataSource;
import org.igv.data.DataTile;
import org.igv.data.NamedScore;
import org.igv.feature.LocusScore;
import org.igv.feature.genome.Genome;
import org.igv.track.DataType;
import org.igv.track.WindowFunction;
import org.igv.ui.panel.FrameManager;
import org.igv.util.collections.LRUCache;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;

/**
 * One track of a figwig file.  Raw data are the track's nonzero bases, merged into runs of one value;
 * zoomed-out views use the tile sums and zoom levels, which hold sums, so only the mean is precomputed.
 */
public class FigWigDataSource extends AbstractDataSource implements DataSource {

    private static final org.igv.logging.Logger log = org.igv.logging.LogManager.getLogger(FigWigDataSource.class);

    private static final Collection<WindowFunction> WINDOW_FUNCTIONS =
            Arrays.asList(WindowFunction.mean, WindowFunction.none);

    private final FigWigFile file;
    private final int track;
    private final Map<String, Integer> chromIds = new HashMap<>();

    public FigWigDataSource(FigWigFile file, int track, Genome genome) {
        super(genome);
        this.file = file;
        this.track = track;
        file.retain(track);
        for (int c = 0; c < file.chromNames.length; c++) {
            String name = file.chromNames[c];
            chromIds.put(name, c);
            if (genome != null) {
                chromIds.putIfAbsent(genome.getCanonicalChrName(name), c);
            }
        }
    }

    public FigWigFile getFile() {
        return file;
    }

    public double getDataMax() {
        return file.statMax[track];
    }

    public double getDataMin() {
        return Math.min(0, file.statMin[track]);
    }

    public DataType getDataType() {
        return DataType.OTHER;
    }

    @Override
    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return WINDOW_FUNCTIONS;
    }

    @Override
    public int getLongestFeature(String chr) {
        return 0;
    }

    @Override
    protected DataTile getRawData(String chr, int start, int end) {
        Integer c = chromIds.get(chr);
        if (c == null) {
            return null;
        }
        try {
            FigWigFile.Runs runs = file.read(track, c, start, end);
            return tile(runs);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    @Override
    protected List<LocusScore> getPrecomputedSummaryScores(String chr, int start, int end, int zoom) {
        if (windowFunction != WindowFunction.mean) {
            return null;
        }
        if (Globals.CHR_ALL.equals(chr)) {
            return wholeGenomeScores();
        }
        Integer c = chromIds.get(chr);
        if (c == null) {
            return Collections.emptyList();
        }
        double scale = file.chromLengths[c] / (Math.pow(2, zoom) * 700);
        if (file.tile > scale || !file.hasSummary(track)) {
            // The finest level (tile sums) is coarser than a bin, or the group keeps none: raw data.
            return null;
        }
        int level = -1;
        for (int k = 0; k < file.levels(); k++) {
            if (file.levelBases(k) <= scale) {
                level = k;
            }
        }
        if (level < 0) {
            return null;
        }
        FigWigFile.Runs runs = file.levelMeans(track, level, c, start, end);
        List<LocusScore> scores = new ArrayList<>(runs.n);
        for (int i = 0; i < runs.n; i++) {
            scores.add(new BasicScore(runs.starts[i], runs.ends[i], (float) runs.values[i]));
        }
        return scores;
    }

    private final Map<WindowFunction, List<LocusScore>> wholeGenome = new HashMap<>();

    /**
     * The whole-genome view, as BBDataSource draws it: one level of summary over every long
     * chromosome, in genome coordinates (kilobases), computed once.  The level is the coarsest whose
     * bins are no larger than the genome over 1,000 nominal pixels.
     */
    private List<LocusScore> wholeGenomeScores() {
        if (genome == null || !Globals.CHR_ALL.equals(genome.getHomeChromosome())) {
            return null;
        }
        synchronized (wholeGenome) {
            List<LocusScore> scores = wholeGenome.get(windowFunction);
            if (scores != null) {
                return scores;
            }
            double scale = genome.getWGLength() / 1000.0;
            int level = -1;
            for (int k = 0; k < file.levels(); k++) {
                if (file.levelBases(k) <= scale) {
                    level = k;
                }
            }
            if (level < 0) {
                return null;
            }
            scores = new ArrayList<>();
            for (String chr : genome.getLongChromosomeNames()) {
                Integer c = chromIds.get(chr);
                if (c == null) {
                    continue;
                }
                FigWigFile.Runs runs = file.levelMeans(track, level, c, 0, (int) file.chromLengths[c]);
                if (runs == null) {
                    break;  // the group keeps no summary: nothing drawn, as reading the genome raw is too slow
                }
                for (int i = 0; i < runs.n; i++) {
                    int a = genome.getGenomeCoordinate(chr, runs.starts[i]);
                    int b = genome.getGenomeCoordinate(chr, runs.ends[i]);
                    scores.add(new BasicScore(a, Math.max(a + 1, b), (float) runs.values[i]));
                }
            }
            scores.sort(Comparator.comparingInt(LocusScore::getStart));
            wholeGenome.put(windowFunction, scores);
            return scores;
        }
    }

    /**
     * The level a zoom would use, or -1 for raw data; for tests and timing.
     */
    public int levelFor(String chr, int zoom) {
        Integer c = chromIds.get(chr);
        double scale = file.chromLengths[c] / (Math.pow(2, zoom) * 700);
        int level = -1;
        for (int k = 0; k < file.levels(); k++) {
            if (file.levelBases(k) <= scale) {
                level = k;
            }
        }
        return level;
    }

    // ------------------------------------------------------------------------------------------
    // IGV's summary of raw data, for the mean and none window functions, computed from the reader's
    // runs without a DataTile or an Accumulator per bin.  It reproduces AbstractDataSource's
    // getSummaryScoresForRange, getSummaryTilesForRange and computeSummaryTile step for step, float
    // arithmetic, score classes and representative values included; other window functions use them.

    private final LRUCache<String, List<LocusScore>> tileCache = new LRUCache<>(10);

    // Loading every tile of a request in one pass first: "http" (the default) for remote files only,
    // where it overlaps round trips; "always"; or "never".
    private static final String PREFETCH = System.getProperty("figwig.prefetch", "http");

    @Override
    public List<LocusScore> getSummaryScoresForRange(String chr, int startLocation, int endLocation, int zoom) {
        if (windowFunction != WindowFunction.mean && windowFunction != WindowFunction.none) {
            return super.getSummaryScoresForRange(chr, startLocation, endLocation, zoom);
        }
        if (windowFunction != WindowFunction.none) {
            List<LocusScore> scores = getPrecomputedSummaryScores(chr, startLocation, endLocation, zoom);
            if (scores != null) {
                return scores;
            }
        }
        List<List<LocusScore>> tiles = summaryTiles(chr, startLocation, endLocation, zoom);
        List<LocusScore> scores = new ArrayList<>(tiles.size() * 700);
        for (List<LocusScore> tile : tiles) {
            scores.addAll(tile);
        }
        return scores;
    }

    /** Requests for raw summaries running now, over every figwig track. */
    private static final AtomicInteger ACTIVE_SUMMARIES = new AtomicInteger();

    private List<List<LocusScore>> summaryTiles(String chr, int startLocation, int endLocation, int zReq) {
        int chrLength = getChrLength(chr);
        if (chrLength == 0) {
            return Collections.emptyList();
        }
        endLocation = Math.min(endLocation, chrLength);
        int adjustedStart = Math.max(0, startLocation);
        int adjustedEnd = Math.min(chrLength, endLocation);
        if (FrameManager.isGeneListMode()) {
            return Collections.singletonList(summarize(chr, startLocation, endLocation, 700));
        }
        int z = zReq;
        int virtualTileCount = (int) Math.pow(2, z);
        double tileWidth = ((double) chrLength) / virtualTileCount;
        int startTile = (int) (adjustedStart / tileWidth);
        int endTile = (int) (Math.min(chrLength, adjustedEnd) / tileWidth) + 1;
        // Over HTTP, load the chunks of every tile not cached in one pass, so that their reads overlap.
        Integer c = chromIds.get(chr);
        if (c != null && (file.isHttp() ? !"never".equals(PREFETCH) : "always".equals(PREFETCH))) {
            int missFirst = -1, missLast = -1;
            for (int t = startTile; t <= endTile; t++) {
                if (tileCache.get(chr + "_" + z + "_" + t + getWindowFunction()) == null) {
                    if (missFirst < 0) missFirst = t;
                    missLast = t;
                }
            }
            if (missFirst >= 0) {
                try {
                    file.prefetch(track, c, (int) (missFirst * tileWidth), Math.min(chrLength, (int) ((missLast + 1) * tileWidth)));
                } catch (IOException e) {
                    throw new UncheckedIOException(e);
                }
            }
        }
        // Large tiles not cached are summarized in parallel, each on its own as IGV does; small ones on
        // this thread, as handing them to the pool costs more than they do.
        if (FigWigPool.DECODE_THREADS > 1 && endTile > startTile && tileWidth >= 8.0 * file.tile) {
            // Summaries large enough to run in parallel are counted, so that those running at once share the threads.
            if (ACTIVE_SUMMARIES.incrementAndGet() * 2 <= FigWigPool.DECODE_THREADS) {
                try {
                    return summaryTilesParallel(chr, chrLength, z, tileWidth, startTile, endTile);
                } finally {
                    ACTIVE_SUMMARIES.decrementAndGet();
                }
            }
            ACTIVE_SUMMARIES.decrementAndGet();
        }
        List<List<LocusScore>> tiles = new ArrayList<>(endTile - startTile + 1);
        for (int t = startTile; t <= endTile; t++) {
            int tileStart = (int) (t * tileWidth);
            int tileEnd = Math.min(chrLength, (int) ((t + 1) * tileWidth));
            String key = chr + "_" + z + "_" + t + getWindowFunction();
            List<LocusScore> tile = tileCache.get(key);
            if (tile == null) {
                tile = summarize(chr, tileStart, tileEnd, 700);
                synchronized (tileCache) {
                    tileCache.put(key, tile);
                }
            }
            tiles.add(tile);
        }
        return tiles;
    }

    private List<List<LocusScore>> summaryTilesParallel(String chr, int chrLength, int z, double tileWidth,
                                                        int startTile, int endTile) {
        int count = endTile - startTile + 1;
        List<LocusScore>[] got = new List[count];
        List<Callable<List<LocusScore>>> tasks = new ArrayList<>();
        List<Integer> missing = new ArrayList<>();
        for (int t = startTile; t <= endTile; t++) {
            got[t - startTile] = tileCache.get(chr + "_" + z + "_" + t + getWindowFunction());
            if (got[t - startTile] == null) {
                final int tileStart = (int) (t * tileWidth);
                final int tileEnd = Math.min(chrLength, (int) ((t + 1) * tileWidth));
                tasks.add(() -> summarize(chr, tileStart, tileEnd, 700));
                missing.add(t);
            }
        }
        if (!tasks.isEmpty()) {
            List<List<LocusScore>> done;
            try {
                done = FigWigPool.parallel(tasks);
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
            for (int i = 0; i < missing.size(); i++) {
                int t = missing.get(i);
                got[t - startTile] = done.get(i);
                synchronized (tileCache) {
                    tileCache.put(chr + "_" + z + "_" + t + getWindowFunction(), done.get(i));
                }
            }
        }
        return Arrays.asList(got);
    }

    /** AbstractDataSource.computeSummaryTile over the reader's runs. */
    private List<LocusScore> summarize(String chr, int startLocation, int endLocation, int nBins) {
        List<LocusScore> out = new ArrayList<>(1000);
        Integer c = chromIds.get(chr);
        if (c == null || nBins <= 0) {
            return out;
        }
        FigWigFile.Runs runs;
        try {
            runs = file.read(track, c, startLocation, endLocation);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
        int size = runs.n;
        if (size == 0) {
            return out;
        }
        int[] starts = runs.starts;
        int[] ends = runs.ends;
        double[] values = runs.values;
        if (windowFunction == WindowFunction.none) {
            for (int i = 0; i < size; i++) {
                int s = starts[i];
                int e = Math.max(s + 1, ends[i]);
                if (e < startLocation) {
                    continue;
                } else if (s >= endLocation) {
                    break;
                }
                out.add(new NamedScore(s, e, (float) values[i], null));
            }
            return out;
        }
        double scale = (double) (endLocation - startLocation) / nBins;
        // One Accumulator's state (mean): its float sum, bases, points and first five values.
        float sum = 0.0f;
        int basesCovered = 0;
        int nPts = 0;
        float[] rep = new float[5];
        int accumulatedStart = -1;
        int accumulatedEnd = -1;
        int lastEndBin = 0;
        for (int i = 0; i < size; i++) {
            int trueEnd = ends[i];
            float v = (float) values[i] * 1.0f;
            if (starts[i] >= endLocation) {
                break;
            } else if (trueEnd <= startLocation || Float.isNaN(v)) {
                continue;
            }
            int s = Math.max(startLocation, starts[i]);
            int e = Math.min(endLocation, trueEnd);
            int endBin = (int) ((e - startLocation) / scale);
            int startBin = (int) ((s - startLocation) / scale);
            if (endBin > lastEndBin || endBin > startBin) {
                if (basesCovered > 0) {
                    out.add(composite(accumulatedStart, accumulatedEnd, sum, basesCovered, nPts, rep));
                    sum = 0.0f;
                    basesCovered = 0;
                    nPts = 0;
                }
            }
            if (endBin > startBin) {
                out.add(new NamedScore(s, e, v, null));
            } else {
                if (basesCovered <= 0) {
                    accumulatedStart = s;
                }
                accumulatedEnd = e;
                int nBases = e - s;
                if (nBases < 1) {
                    nBases = 1;
                }
                if (nPts < 5) {
                    rep[nPts] = v;
                }
                sum += nBases * v;
                nPts++;
                basesCovered += nBases;
            }
            lastEndBin = endBin;
        }
        if (basesCovered > 0) {
            out.add(composite(accumulatedStart, accumulatedEnd, sum, basesCovered, nPts, rep));
        }
        return out;
    }

    /** AbstractDataSource.getCompositeScore for a mean Accumulator in this state. */
    private LocusScore composite(int start, int end, float sum, int basesCovered, int nPts, float[] rep) {
        if (nPts == 1) {
            return new NamedScore(start, end, rep[0], null);
        }
        float value = Float.isNaN(sum) ? Float.NaN : sum / basesCovered;
        float[] data = new float[5];
        System.arraycopy(rep, 0, data, 0, Math.min(nPts, 5));
        return new CompositeScore(start, end, value, data, new String[5], windowFunction);
    }

    private static DataTile tile(FigWigFile.Runs runs) {
        int[] starts = Arrays.copyOf(runs.starts, runs.n);
        int[] ends = Arrays.copyOf(runs.ends, runs.n);
        float[] values = new float[runs.n];
        for (int i = 0; i < runs.n; i++) {
            values[i] = (float) runs.values[i];
        }
        return new DataTile(starts, ends, values, null);
    }

    private boolean disposed = false;

    /** IGV disposes a track's source when the track is removed; the last track of a file closes it. */
    @Override
    public synchronized void dispose() {
        if (!disposed) {
            disposed = true;
            try {
                file.release(track);
            } catch (IOException e) {
                log.error("Error closing " + file.path, e);
            }
        }
    }
}
