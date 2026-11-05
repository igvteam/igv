package org.igv.seg;

import org.igv.feature.LocusScore;
import org.igv.track.DataType;
import org.igv.track.WindowFunction;

import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public interface SegmentedDataSource {

    DataType getType();

    double getDataMin();

    double getDataMax();

    List<String> getSampleNames();

    List<LocusScore> getFeatures(String sample, String chr);

    /**
     * Scores for every sample over a region, keyed by sample name.  In-memory sources ignore the range and zoom
     * and return their per-chromosome lists; summarizing sources fetch the region once for all samples.
     */
    default Map<String, List<LocusScore>> getFeatures(String chr, int start, int end, int zoom) {
        Map<String, List<LocusScore>> features = new LinkedHashMap<>();
        for (String sample : getSampleNames()) {
            features.put(sample, getFeatures(sample, chr));
        }
        return features;
    }

    default void setWindowFunction(WindowFunction windowFunction) {
    }

    default WindowFunction getWindowFunction() {
        return null;
    }

    default Collection<WindowFunction> getAvailableWindowFunctions() {
        return Collections.emptyList();
    }

    default boolean isLogNormalized() {
        return true;
    }

    default void dispose() {
    }
}
