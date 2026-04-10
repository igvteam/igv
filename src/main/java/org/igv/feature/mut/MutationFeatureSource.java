package org.igv.feature.mut;

import htsjdk.tribble.readers.AsciiLineReader;
import org.igv.exceptions.DataLoadException;
import org.igv.feature.LocusScore;
import org.igv.feature.genome.Genome;
import org.igv.feature.tribble.MUTCodec;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.seg.SegmentedDataSource;
import org.igv.track.DataType;
import org.igv.track.FeatureCollectionSource;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.util.ParsingUtils;
import org.igv.util.ResourceLocator;

import java.io.IOException;
import java.util.*;

/**
 * @author jrobinso
 *         Date: 4/9/13
 *         Time: 8:45 AM
 */
public class MutationFeatureSource implements SegmentedDataSource {

    private static final Logger log = LogManager.getLogger(MutationFeatureSource.class);

    private final ResourceLocator locator;
    private final Genome genome;
    // Map of sampleId -> FeatureCollectionSource<Mutation>
    private Map<String, FeatureCollectionSource<Mutation>> featureCollectionSourceMap;
    List<String> samples;

    private double dataMin = Double.POSITIVE_INFINITY;
    private double dataMax = Double.NEGATIVE_INFINITY;

    public MutationFeatureSource(ResourceLocator locator, Genome genome) {
        this.locator = locator;
        this.genome = genome;
        loadMutations();
    }


    @Override
    public List<String> getSampleNames() {
        return samples;
    }

    @Override
    public List<LocusScore> getFeatures(String sample, String chr) {
        FeatureCollectionSource<Mutation> source = featureCollectionSourceMap.get(sample);
        if (source == null) return List.of();
        List<LocusScore> result = new ArrayList<>();
        Iterator<Mutation> it = source.getFeatures(chr, 0, Integer.MAX_VALUE);
        while (it.hasNext()) {
            result.add(it.next());
        }
        return result;
    }

    @Override
    public LocusScore getFeatureAt(String sample, String chr, double position, ReferenceFrame frame) {
        FeatureCollectionSource<Mutation> source = featureCollectionSourceMap.get(sample);
        if (source == null) return null;
        // Add a 3 pixel tolerance
        double locScale = frame.getScale();
        int tolerance = (int) (3 * locScale);
        Iterator<Mutation> it = source.getFeatures(chr, (int) (position - tolerance), (int) (position + tolerance));
        Mutation closest = null;
        double minDist = Double.POSITIVE_INFINITY;
        while (it.hasNext()) {
            Mutation mut = it.next();
            double dist = Math.abs(mut.getStart() - position);
            if (dist < minDist) {
                minDist = dist;
                closest = mut;
            }
        }
        return closest;
    }

    /**
     * Load mutations and group by sampleId as they are read.
     */
    private void loadMutations() {
        AsciiLineReader reader = null;
        samples = new ArrayList<>();
        Set<String> sampleSet = new HashSet<>();
        featureCollectionSourceMap = new LinkedHashMap<>();
        Map<String, List<Mutation>> mutationMap = new LinkedHashMap<>();
        try {
            MUTCodec codec = new MUTCodec(locator.getPath(), genome);
            reader = ParsingUtils.openAsciiReader(locator);
            String nextLine;
            // Skip header - handled in codec
            while ((nextLine = reader.readLine()) != null) {
                if (nextLine.startsWith("#")) {
                    continue;
                } else {
                    break;
                }
            }
            // Read first data line if not null
            if (nextLine != null) {
                Mutation mut = codec.decode(nextLine);
                if (mut != null) {
                    String sampleId = mut.getSampleId();
                    if (sampleSet.add(sampleId)) {
                        samples.add(sampleId);
                    }
                    mutationMap.computeIfAbsent(sampleId, k -> new ArrayList<>()).add(mut);
                    double score = mut.getScore();
                    if (!Double.isNaN(score)) {
                        dataMin = Math.min(dataMin, score);
                        dataMax = Math.max(dataMax, score);
                    }
                }
            }
            while ((nextLine = reader.readLine()) != null) {
                Mutation mut = codec.decode(nextLine);
                if (mut != null) {
                    String sampleId = mut.getSampleId();
                    if (sampleSet.add(sampleId)) {
                        samples.add(sampleId);
                    }
                    mutationMap.computeIfAbsent(sampleId, k -> new ArrayList<>()).add(mut);
                    double score = mut.getScore();
                    if (!Double.isNaN(score)) {
                        dataMin = Math.min(dataMin, score);
                        dataMax = Math.max(dataMax, score);
                    }
                }
            }
            // Now create FeatureCollectionSource for each sample
            for (Map.Entry<String, List<Mutation>> entry : mutationMap.entrySet()) {
                featureCollectionSourceMap.put(entry.getKey(), new FeatureCollectionSource<>(entry.getValue(), genome));
            }
        } catch (IOException e) {
            log.error("Error loading mutation file", e);
            throw new DataLoadException("IO Exception: " + e, locator.getPath());
        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (Exception ignored) {
                }
            }
        }
    }

    @Override
    public DataType getType() {
        return DataType.MUTATION;
    }

    @Override
    public double getDataMin() {
        return dataMin == Double.POSITIVE_INFINITY ? 0 : dataMin;
    }

    @Override
    public double getDataMax() {
        return dataMax == Double.NEGATIVE_INFINITY ? 0 : dataMax;
    }

}
