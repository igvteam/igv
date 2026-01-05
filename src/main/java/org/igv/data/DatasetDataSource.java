package org.igv.data;


import org.igv.logging.*;
import org.igv.Globals;
import org.igv.feature.Chromosome;
import org.igv.feature.LocusScore;
import org.igv.feature.genome.Genome;
import org.igv.track.TrackType;
import org.igv.track.WindowFunction;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 */
public class DatasetDataSource extends AbstractDataSource {

    private Logger log = LogManager.getLogger(DatasetDataSource.class);

    String trackId;
    Dataset dataset;
    GenomeSummaryData genomeSummaryData;

    /**
     * @param trackId
     * @param dataset
     * @param genome
     */
    public DatasetDataSource(String trackId, Dataset dataset, Genome genome) {
        super(genome);
        this.trackId = trackId;
        this.dataset = dataset;

        if (genome != null && genome.getLongChromosomeNames() != null && genome.getLongChromosomeNames().size() > 0) {
            if (dataset instanceof IGVDataset) {
                genomeSummaryData = ((IGVDataset) dataset).getGenomeSummary();
            } else {
                genomeSummaryData = new GenomeSummaryData(genome, new String[]{trackId});
                for (String chrName : genome.getLongChromosomeNames()) {
                    Chromosome chr = genome.getChromosome(chrName);
                    int[] startLocations = dataset.getStartLocations(chr.getName());
                    if (!chr.getName().equals(Globals.CHR_ALL) && (startLocations != null) && (startLocations.length > 0)) {
                        Map<String, float[]> dMap = new HashMap<String, float[]>();
                        dMap.put(trackId, dataset.getData(trackId, chr.getName()));
                        genomeSummaryData.addData(chr.getName(), startLocations, dMap);
                    }
                }
            }

        }
    }


    @Override
    protected DataTile getRawData(String chr, int startLocation, int endLocation) {

        if (chr.equals(Globals.CHR_ALL) && genomeSummaryData != null && windowFunction != WindowFunction.none) {
            int[] startLocs = genomeSummaryData.getLocations();
            int[] endLocs = null;
            float[] data = genomeSummaryData.getData(trackId);
            return new DataTile(startLocs, endLocs, data, null);
        }
        if (chr.equals(Globals.CHR_ALL)) {
            return getWGRawData();
        } else {
            int[] startLocs = dataset.getStartLocations(chr);
            int[] endLocs = dataset.getEndLocations(chr);
            float[] data = dataset.getData(trackId, chr);
            String[] features = dataset.getFeatureNames(chr);

            if (startLocs == null || (data == null) || data.length == 0) {
                return null;
            }

            return new DataTile(startLocs, endLocs, data, features);
        }
    }

    /**
     * Create a tile for WG data.  This is used when genome summary data is not useful => window function == none
     *
     * @return
     */
    private DataTile getWGRawData() {

        int size = 0;
        for (String chr : genome.getChromosomeNames()) {
            int[] s = dataset.getStartLocations(chr);
            int[] e = dataset.getEndLocations(chr);
            float[] d = dataset.getData(trackId, chr);
            if (s != null && d != null) size += s.length;
        }
        if (size == 0) return null;

        int[] startLocs = new int[size];
        int[] endLocs = new int[size];
        float[] data = new float[size];
        String[] features = new String[size];

        int i = 0;
        for (String chr : genome.getChromosomeNames()) {
            int[] s = dataset.getStartLocations(chr);
            int[] e = dataset.getEndLocations(chr);
            float[] d = dataset.getData(trackId, chr);
            String[] f = dataset.getFeatureNames(chr);

            if (s != null && d != null) {
                int l = s.length;
                for (int j = 0; j < l; j++) {

                    startLocs[i] = genome.getGenomeCoordinate(chr, s[j]);
                    endLocs[i] = e == null ? startLocs[i] + 1 : genome.getGenomeCoordinate(chr, e[j]);
                    data[i] = d[j];
                    if (f != null) features[i] = f[j];

                    i++;
                }
            }

        }
        return new DataTile(startLocs, endLocs, data, features);

    }

    @Override
    protected List<LocusScore> getPrecomputedSummaryScores(String chr, int startLocation, int endLocation, int zoom) {
        return null;
    }


    public TrackType getTrackType() {
        try {
            return dataset.getType();
        } catch (Exception exception) {
            return TrackType.OTHER;
        }
    }


    @Override
    public boolean isLogNormalized() {
        return dataset.isLogNormalized();
    }


    public double getDataMax() {
        return dataset.getDataMax();
    }


    public double getDataMin() {
        return dataset.getDataMin();
    }

    @Override
    public int getLongestFeature(String chr) {
        return dataset.getLongestFeature(chr);
    }
}
