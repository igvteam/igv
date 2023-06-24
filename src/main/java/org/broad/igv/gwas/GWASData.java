package org.broad.igv.gwas;

import java.util.*;

/**
 * Simple map-like class that contains GWAS pvalue data organized by chromosome, along with min and max values*
 */

public class GWASData {

    double minValue;
    double maxValue;
    Map<String, List<GWASFeature>> features;

    public GWASData() {
        this.features = new HashMap<>();
        minValue = Double.MAX_VALUE;
        maxValue = - Double.MAX_VALUE;
    }

    void addFeature(GWASFeature f) {
        List<GWASFeature> featureList = features.get(f.chr);
        if (featureList == null) {
            featureList = new ArrayList<>();
            features.put(f.chr, featureList);
        }
        featureList.add(f);
        minValue = Math.min(minValue, f.value);
        maxValue = Math.max(maxValue, f.value);
    }

    Set<String> keySet() {
        return features.keySet();
    }
    Collection<List<GWASFeature>> values() {
        return features.values();
    }

    List<GWASFeature> get(String chr) {
        return features.get(chr);
    }

    boolean isEmpty() {
        return features.isEmpty();
    }

    public double getMinValue() {
        return minValue;
    }

    public double getMaxValue() {
        return maxValue;
    }

    void finish() {
        // Sort features by position
        for (List<GWASFeature> featureList : features.values()) {
            featureList.sort(Comparator.comparingInt(o -> o.position));
        }
    }
}
