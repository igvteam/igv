package org.igv.data;

import org.igv.Globals;
import org.igv.feature.Chromosome;
import org.igv.feature.LocusScore;
import org.igv.feature.genome.Genome;
import org.igv.seg.SegmentedDataSource;
import org.igv.tdf.Accumulator;
import org.igv.track.WindowFunction;
import org.igv.util.collections.LRUCache;

import java.util.*;

/**
 * A numeric data source for several samples sharing one set of genomic positions (a locus-by-sample matrix).  Data
 * for a region is fetched and summarized once for all samples, then split into per-sample score lists.  Subclasses
 * tile the genome as the single-sample source they mirror does; tiles are cached here, each holding every sample.
 */
public abstract class MultiSampleDataSource implements SegmentedDataSource {

    private static final int BINS_PER_TILE = 700;

    protected final Genome genome;
    protected final List<String> sampleNames;
    protected WindowFunction windowFunction = WindowFunction.mean;
    private final LRUCache<String, Map<String, List<LocusScore>>> tileCache = new LRUCache<>(10);

    protected MultiSampleDataSource(Genome genome, List<String> sampleNames) {
        this.genome = genome;
        this.sampleNames = sampleNames;
    }

    /**
     * Compute the scores of every sample for the tile [start, end) at the given zoom.
     */
    protected abstract Map<String, List<LocusScore>> computeTile(String chr, int start, int end, int zoom);

    @Override
    public abstract Map<String, List<LocusScore>> getFeatures(String chr, int start, int end, int zoom);

    @Override
    public List<String> getSampleNames() {
        return sampleNames;
    }

    protected int getChrLength(String chr) {
        if (chr.equals(Globals.CHR_ALL)) {
            return (int) (genome.getWGLength() / 1000);
        } else {
            Chromosome c = genome.getChromosome(chr);
            return c == null ? 0 : c.getLength();
        }
    }

    protected Map<String, List<LocusScore>> getCachedTile(String chr, int zoom, int tileNumber, int tileStart, int tileEnd) {
        String key = chr + "_" + zoom + "_" + tileNumber + "_" + windowFunction;
        Map<String, List<LocusScore>> tile = tileCache.get(key);
        if (tile == null) {
            tile = computeTile(chr, tileStart, tileEnd, zoom);
            synchronized (tileCache) {
                tileCache.put(key, tile);
            }
        }
        return tile;
    }

    protected static void addAll(Map<String, List<LocusScore>> target, Map<String, List<LocusScore>> tile) {
        for (Map.Entry<String, List<LocusScore>> entry : tile.entrySet()) {
            List<LocusScore> list = target.get(entry.getKey());
            if (list != null) {
                list.addAll(entry.getValue());
            }
        }
    }

    protected Map<String, List<LocusScore>> emptyScores() {
        Map<String, List<LocusScore>> scores = new LinkedHashMap<>();
        for (String sample : sampleNames) {
            scores.put(sample, new ArrayList<>());
        }
        return scores;
    }

    /**
     * The lists of a scores map as an array in sampleNames order.
     */
    protected List<LocusScore>[] sampleLists(Map<String, List<LocusScore>> scores) {
        List<LocusScore>[] lists = new List[sampleNames.size()];
        for (int n = 0; n < lists.length; n++) {
            lists[n] = scores.get(sampleNames.get(n));
        }
        return lists;
    }

    @Override
    public List<LocusScore> getFeatures(String sample, String chr) {
        return getFeatures(chr, 0, getChrLength(chr), 0).get(sample);
    }

    @Override
    public void setWindowFunction(WindowFunction windowFunction) {
        this.windowFunction = windowFunction;
        tileCache.clear();
    }

    @Override
    public WindowFunction getWindowFunction() {
        return windowFunction;
    }

    /**
     * Bin the values of every sample over [startLocation, endLocation) into BINS_PER_TILE bins in one pass.  The bin
     * arithmetic for each feature is done once; each sample keeps its own accumulator state so the per-sample result
     * equals what AbstractDataSource.computeSummaryTile produces for that sample alone.
     *
     * @param values one array per sample, in sampleNames order, aligned with starts.  A null entry yields no scores.
     */
    protected Map<String, List<LocusScore>> summarize(int[] starts, int[] ends, String[] names, float[][] values,
                                                      int startLocation, int endLocation) {

        Map<String, List<LocusScore>> result = emptyScores();
        if (starts == null || starts.length == 0) {
            return result;
        }
        int nSamples = sampleNames.size();
        List<LocusScore>[] lists = sampleLists(result);

        if (windowFunction == WindowFunction.none) {

            for (int i = 0; i < starts.length; i++) {
                int s = starts[i];
                int e = ends == null ? s + 1 : Math.max(s + 1, ends[i]);
                if (e < startLocation) {
                    continue;
                } else if (s >= endLocation) {
                    break;
                }
                String probeName = names == null ? null : names[i];
                for (int n = 0; n < nSamples; n++) {
                    if (values[n] != null) {
                        lists[n].add(new NamedScore(s, e, values[n][i], probeName));
                    }
                }
            }

        } else {

            double scale = (double) (endLocation - startLocation) / BINS_PER_TILE;
            Accumulator[] accumulators = new Accumulator[nSamples];
            int[] accumulatedStart = new int[nSamples];
            int[] accumulatedEnd = new int[nSamples];
            int[] lastEndBin = new int[nSamples];

            for (int i = 0; i < starts.length; i++) {

                int trueEnd = ends == null ? starts[i] + 1 : ends[i];
                if (starts[i] >= endLocation) {
                    break;  // We're beyond the end of the requested interval
                } else if (trueEnd <= startLocation) {
                    continue;
                }

                // Bound feature at interval, other "piece" will be in another tile.
                int s = Math.max(startLocation, starts[i]);
                int e = Math.min(endLocation, trueEnd);
                String probeName = names == null ? null : names[i];

                // Compute bin numbers, relative to start of this tile
                int endBin = (int) ((e - startLocation) / scale);
                int startBin = (int) ((s - startLocation) / scale);

                for (int n = 0; n < nSamples; n++) {
                    if (values[n] == null) {
                        continue;
                    }
                    float v = values[n][i];
                    if (Float.isNaN(v)) {
                        continue;
                    }

                    // If this feature spans multiple bins, or extends beyond last end bin, record
                    if (endBin > lastEndBin[n] || endBin > startBin) {
                        if (accumulators[n] != null && accumulators[n].hasData()) {
                            lists[n].add(getCompositeScore(accumulators[n], accumulatedStart[n], accumulatedEnd[n]));
                            accumulators[n] = null;
                        }
                    }

                    if (endBin > startBin) {
                        lists[n].add(new NamedScore(s, e, v, probeName));
                    } else {
                        if (accumulators[n] == null) {
                            accumulators[n] = new Accumulator(windowFunction, 5);
                        }
                        if (!accumulators[n].hasData()) accumulatedStart[n] = s;
                        accumulatedEnd[n] = e;
                        accumulators[n].add(e - s, v, probeName);
                    }

                    lastEndBin[n] = endBin;
                }
            }

            // Cleanup
            for (int n = 0; n < nSamples; n++) {
                if (accumulators[n] != null && accumulators[n].hasData()) {
                    lists[n].add(getCompositeScore(accumulators[n], accumulatedStart[n], accumulatedEnd[n]));
                }
            }
        }

        return result;
    }

    private LocusScore getCompositeScore(Accumulator accumulator, int accumulatedStart, int accumulatedEnd) {
        LocusScore ls;
        if (accumulator.getNpts() == 1) {
            ls = new NamedScore(accumulatedStart, accumulatedEnd, accumulator.getRepData()[0], accumulator.getRepProbes()[0]);
        } else {
            float value = accumulator.getValue();
            ls = new CompositeScore(accumulatedStart, accumulatedEnd, value, accumulator.getRepData(),
                    accumulator.getRepProbes(), windowFunction);
        }
        return ls;
    }
}
