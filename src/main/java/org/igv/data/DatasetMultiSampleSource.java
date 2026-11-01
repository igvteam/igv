package org.igv.data;

import org.igv.Globals;
import org.igv.feature.Chromosome;
import org.igv.feature.LocusScore;
import org.igv.feature.genome.Genome;
import org.igv.track.DataType;
import org.igv.track.WindowFunction;

import java.util.*;

/**
 * Multi-sample source backed by an in-memory Dataset (.igv, .cn, .gct, ...).  The dataset's positions for a
 * chromosome are shared by all samples, so a tile is summarized in one pass over them.
 */
public class DatasetMultiSampleSource extends MultiSampleDataSource {

    private final Dataset dataset;
    private GenomeSummaryData genomeSummaryData;

    public DatasetMultiSampleSource(Dataset dataset, Genome genome) {
        super(genome, Arrays.asList(dataset.getSampleNames()));
        this.dataset = dataset;

        if (genome != null && genome.getLongChromosomeNames() != null && genome.getLongChromosomeNames().size() > 0) {
            if (dataset instanceof IGVDataset) {
                genomeSummaryData = ((IGVDataset) dataset).getGenomeSummary();
            } else {
                genomeSummaryData = new GenomeSummaryData(genome, dataset.getSampleNames());
                for (String chrName : genome.getLongChromosomeNames()) {
                    Chromosome chr = genome.getChromosome(chrName);
                    int[] startLocations = dataset.getStartLocations(chr.getName());
                    if (!chr.getName().equals(Globals.CHR_ALL) && (startLocations != null) && (startLocations.length > 0)) {
                        Map<String, float[]> dMap = new HashMap<>();
                        for (String sample : sampleNames) {
                            dMap.put(sample, dataset.getData(sample, chr.getName()));
                        }
                        genomeSummaryData.addData(chr.getName(), startLocations, dMap);
                    }
                }
            }
        }
    }

    /**
     * Tiling as in AbstractDataSource.getSummaryTilesForRange: 2^zoom tiles per chromosome.
     */
    @Override
    public Map<String, List<LocusScore>> getFeatures(String chr, int start, int end, int zoom) {

        Map<String, List<LocusScore>> scores = emptyScores();
        int chrLength = getChrLength(chr);
        if (chrLength == 0) {
            return scores;
        }
        int adjustedStart = Math.max(0, start);
        int adjustedEnd = Math.min(chrLength, end);

        if (chr.equals(Globals.CHR_ALL)) {
            addAll(scores, getCachedTile(chr, 0, 0, 0, chrLength));
        } else {
            int virtualTileCount = (int) Math.pow(2, zoom);
            double tileWidth = ((double) chrLength) / virtualTileCount;
            int startTile = (int) (adjustedStart / tileWidth);
            int endTile = (int) (adjustedEnd / tileWidth) + 1;
            for (int t = startTile; t <= endTile; t++) {
                int tileStart = (int) (t * tileWidth);
                int tileEnd = Math.min(chrLength, (int) ((t + 1) * tileWidth));
                if (tileStart >= chrLength) {
                    break;
                }
                addAll(scores, getCachedTile(chr, zoom, t, tileStart, tileEnd));
            }
        }
        return scores;
    }

    @Override
    protected Map<String, List<LocusScore>> computeTile(String chr, int start, int end, int zoom) {

        if (chr.equals(Globals.CHR_ALL)) {
            if (genomeSummaryData != null && windowFunction != WindowFunction.none) {
                int[] locations = genomeSummaryData.getLocations();   // Must precede getData, which does not build the arrays
                float[][] values = new float[sampleNames.size()][];
                for (int n = 0; n < values.length; n++) {
                    values[n] = genomeSummaryData.getData(sampleNames.get(n));
                }
                return summarize(locations, null, null, values, start, end);
            }
            return computeWholeGenomeRaw(start, end);
        }

        int[] starts = dataset.getStartLocations(chr);
        if (starts == null || starts.length == 0) {
            return emptyScores();
        }
        return summarize(starts, dataset.getEndLocations(chr), dataset.getFeatureNames(chr), sampleValues(chr), start, end);
    }

    /**
     * The dataset's values for a chromosome, one array per sample in sampleNames order.
     */
    private float[][] sampleValues(String chr) {
        float[][] values = new float[sampleNames.size()][];
        for (int n = 0; n < values.length; n++) {
            values[n] = dataset.getData(sampleNames.get(n), chr);
        }
        return values;
    }

    /**
     * Whole genome data from the raw values of every chromosome, used when there is no genome summary or the
     * window function is "none".
     */
    private Map<String, List<LocusScore>> computeWholeGenomeRaw(int start, int end) {

        int size = 0;
        for (String chr : genome.getChromosomeNames()) {
            int[] s = dataset.getStartLocations(chr);
            if (s != null) size += s.length;
        }
        if (size == 0) return emptyScores();

        int nSamples = sampleNames.size();
        int[] startLocs = new int[size];
        int[] endLocs = new int[size];
        String[] features = new String[size];
        float[][] values = new float[nSamples][size];

        int i = 0;
        for (String chr : genome.getChromosomeNames()) {
            int[] s = dataset.getStartLocations(chr);
            if (s == null) continue;
            int[] e = dataset.getEndLocations(chr);
            String[] f = dataset.getFeatureNames(chr);
            float[][] d = sampleValues(chr);
            for (int j = 0; j < s.length; j++) {
                startLocs[i] = genome.getGenomeCoordinate(chr, s[j]);
                endLocs[i] = e == null ? startLocs[i] + 1 : genome.getGenomeCoordinate(chr, e[j]);
                if (f != null) features[i] = f[j];
                for (int n = 0; n < nSamples; n++) {
                    values[n][i] = d[n] == null ? Float.NaN : d[n][j];
                }
                i++;
            }
        }
        return summarize(startLocs, endLocs, features, values, start, end);
    }

    @Override
    public DataType getType() {
        return dataset.getType();
    }

    @Override
    public double getDataMin() {
        return dataset.getDataMin();
    }

    @Override
    public double getDataMax() {
        return dataset.getDataMax();
    }

    @Override
    public boolean isLogNormalized() {
        return dataset.isLogNormalized();
    }

    @Override
    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return AbstractDataSource.ORDERED_WINDOW_FUNCTIONS;
    }
}
