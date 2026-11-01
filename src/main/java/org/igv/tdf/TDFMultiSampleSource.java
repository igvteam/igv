package org.igv.tdf;

import org.igv.Globals;
import org.igv.data.BasicScore;
import org.igv.data.MultiSampleDataSource;
import org.igv.feature.Chromosome;
import org.igv.feature.LocusScore;
import org.igv.feature.genome.ChromAliasManager;
import org.igv.feature.genome.Genome;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.track.DataType;
import org.igv.track.WindowFunction;

import java.util.*;

/**
 * Multi-sample source backed by a multi-track TDF file.  TDF tiles hold every track, so a tile is read and walked
 * once for all samples rather than once per track as TDFDataSource does.
 */
public class TDFMultiSampleSource extends MultiSampleDataSource {

    private static Logger log = LogManager.getLogger(TDFMultiSampleSource.class);

    private final TDFReader reader;
    private int maxPrecomputedZoom = 6;
    private final Set<String> sequenceNames;
    private final ChromAliasManager chromAliasManager;
    private final List<WindowFunction> availableFunctions;

    public TDFMultiSampleSource(TDFReader reader, Genome genome) {
        super(genome, Arrays.asList(reader.getSampleNames()));
        this.reader = reader;

        sequenceNames = new HashSet<>(reader.getChromosomeNames());
        chromAliasManager = new ChromAliasManager(sequenceNames, genome);
        availableFunctions = new ArrayList<>(reader.getWindowFunctions());
        availableFunctions.add(WindowFunction.none);   // Always available => raw data

        TDFGroup rootGroup = reader.getGroup("/");
        try {
            maxPrecomputedZoom = Integer.parseInt(rootGroup.getAttribute("maxZoom"));
        } catch (Exception e) {
            log.error("Error reading attribute 'maxZoom'", e);
        }
    }

    /**
     * Tiling as in TDFDataSource.getSummaryScoresForRange: the file's tile grid, and only scores overlapping the
     * requested range.
     */
    @Override
    public Map<String, List<LocusScore>> getFeatures(String chr, int start, int end, int zoom) {

        Map<String, List<LocusScore>> scores = emptyScores();
        int chrLength = getChrLength(chr);
        if (chrLength == 0) {
            return scores;
        }
        int startLocation = Math.max(0, start);
        int endLocation = Math.min(chrLength, end);

        int tileWidth = chr.equals(Globals.CHR_ALL) ?
                (int) Math.ceil(genome.getWGLength() / 1000.0) :
                chrLength / ((int) Math.pow(2.0, zoom));
        if (tileWidth <= 0) {
            return scores;
        }
        int startTile = startLocation / tileWidth;
        int endTile = (endLocation - 1) / tileWidth;
        for (int t = startTile; t <= endTile; t++) {
            Map<String, List<LocusScore>> tile = getCachedTile(chr, zoom, t, t * tileWidth, (t + 1) * tileWidth);
            for (Map.Entry<String, List<LocusScore>> entry : tile.entrySet()) {
                List<LocusScore> list = scores.get(entry.getKey());
                for (LocusScore s : entry.getValue()) {
                    if (s.getStart() > endLocation) {
                        break;
                    }
                    if (s.getEnd() >= startLocation) {
                        list.add(s);
                    }
                }
            }
        }
        return scores;
    }

    @Override
    protected Map<String, List<LocusScore>> computeTile(String chr, int start, int end, int zoom) {

        String querySeq = sequenceNames.contains(chr) ? chr : chromAliasManager.getAliasName(chr);
        if (querySeq == null) {
            return emptyScores();
        }

        if (zoom <= maxPrecomputedZoom && windowFunction != WindowFunction.none) {
            List<TDFTile> tiles = null;
            if (querySeq.equals(Globals.CHR_ALL) && !TDFDataSource.isChrOrderValid(reader, genome)) {
                tiles = Arrays.asList(reader.getWholeGenomeTile(genome, windowFunction));
            } else {
                TDFDataset ds = reader.getDataset(querySeq, zoom, windowFunction);
                if (ds != null) {
                    tiles = ds.getTiles(start, end);
                }
            }
            return precomputedScores(tiles);
        } else if (querySeq.equals(Globals.CHR_ALL)) {
            return wholeGenomeRawScores();
        } else {
            return rawScores(querySeq, start, end);
        }
    }

    /**
     * Scores from precomputed (already binned) tiles: one BasicScore per track per bin.
     */
    private Map<String, List<LocusScore>> precomputedScores(List<TDFTile> tiles) {
        Map<String, List<LocusScore>> result = emptyScores();
        if (tiles == null) {
            return result;
        }
        List<LocusScore>[] lists = sampleLists(result);
        for (TDFTile tile : tiles) {
            for (int i = 0; i < tile.getSize(); i++) {
                int s = tile.getStartPosition(i);
                int e = tile.getEndPosition(i);
                for (int n = 0; n < lists.length; n++) {
                    float v = tile.getValue(n, i);
                    if (!Float.isNaN(v)) {
                        lists[n].add(new BasicScore(s, e, v));
                    }
                }
            }
        }
        return result;
    }

    /**
     * Scores from the raw data of a chromosome, binned with the window function (or not, for "none").
     */
    private Map<String, List<LocusScore>> rawScores(String querySeq, int startLocation, int endLocation) {

        Map<String, List<LocusScore>> result = emptyScores();
        TDFDataset rawDataset = reader.getDataset("/" + querySeq + "/raw");
        if (rawDataset == null) {
            return result;
        }
        for (TDFTile rawTile : rawDataset.getTiles(startLocation, endLocation)) {
            if (rawTile == null || rawTile.getSize() == 0) {
                continue;
            }
            float[][] values = new float[sampleNames.size()][];
            for (int n = 0; n < values.length; n++) {
                values[n] = rawTile.getData(n);
            }
            addAll(result, summarize(rawTile.getStart(), rawTile.getEnd(), rawTile.getNames(), values, startLocation, endLocation));
        }
        return result;
    }

    private Map<String, List<LocusScore>> wholeGenomeRawScores() {

        Map<String, List<LocusScore>> result = emptyScores();
        List<LocusScore>[] lists = sampleLists(result);
        for (String chr : genome.getChromosomeNames()) {
            Chromosome c = genome.getChromosome(chr);
            TDFDataset rawDataset = reader.getDataset("/" + chr + "/raw");
            if (rawDataset == null) {
                continue;
            }
            for (TDFTile rawTile : rawDataset.getTiles(0, c.getLength())) {
                if (rawTile == null || rawTile.getSize() == 0) {
                    continue;
                }
                for (int i = 0; i < rawTile.getSize(); i++) {
                    int s = genome.getGenomeCoordinate(chr, rawTile.getStartPosition(i));
                    int e = genome.getGenomeCoordinate(chr, Math.max(rawTile.getStartPosition(i), rawTile.getEndPosition(i) - 1));
                    for (int n = 0; n < lists.length; n++) {
                        lists[n].add(new BasicScore(s, e, rawTile.getValue(n, i)));
                    }
                }
            }
        }
        return result;
    }

    @Override
    public DataType getType() {
        return reader.getTrackType();
    }

    @Override
    public double getDataMin() {
        return reader.getLowerLimit();
    }

    @Override
    public double getDataMax() {
        return reader.getUpperLimit();
    }

    @Override
    public boolean isLogNormalized() {
        return getDataMin() < 0;
    }

    @Override
    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return availableFunctions;
    }

    @Override
    public void dispose() {
        reader.close();
    }
}
