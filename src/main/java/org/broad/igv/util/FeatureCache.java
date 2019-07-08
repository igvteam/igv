package org.broad.igv.util;

import htsjdk.tribble.Feature;
import org.broad.igv.feature.genome.ChromosomeNameComparator;

import java.util.*;

public class FeatureCache {

    int DEFAULT_BATCH_SIZE = 100;

    int batchSize;
    Map<String, IntervalTree<List<Feature>>> featureMap;

    public FeatureCache(List<Feature> features) {

        init(features, DEFAULT_BATCH_SIZE);

    }

    public FeatureCache(List<Feature> features, int batchSize) {
        init(features, batchSize);
    }

    public List<Feature> getFeatures(String chr, int start, int end) {
        List<Feature> features = new ArrayList<>();
        IntervalTree<List<Feature>> tree = featureMap.get(chr);
        if(tree != null) {
            List<Interval<List<Feature>>> intervals = tree.findOverlapping(start, end);
            for (Interval<List<Feature>> interval : intervals) {
                List<Feature> intervalFeatures = interval.getValue();
                for (Feature f : intervalFeatures) {
                    if (f.getEnd() >= start && f.getStart() <= end) {
                        features.add(f);
                    }
                }
            }
        }
        return features;
    }

    private void init(List<Feature> features, int batchSize) {
        String lastChr = null;

        List<Feature> currentFeatureList = new ArrayList<>();
        int currentMin = Integer.MAX_VALUE;
        int currentMax = 0;

        featureMap = new HashMap<>();

        features.sort(getPositionComparator());


        for (Feature f : features) {

            final String chr = f.getChr();
            final int start = f.getStart();
            final int end = f.getEnd();

            if (lastChr == null) {
                currentMin = start;
                currentMax = end;
                currentFeatureList.add(f);
                IntervalTree<List<Feature>> tree = new IntervalTree<>();
                featureMap.put(chr, tree);
                lastChr = chr;
            } else {

                if (!chr.equals(lastChr)) {

                    // New tree
                    IntervalTree<List<Feature>> tree = featureMap.get(lastChr);
                    tree.insert(new Interval(currentMin, currentMax, currentFeatureList));

                    tree = new IntervalTree<>();
                    featureMap.put(chr, tree);
                    lastChr = chr;

                    currentFeatureList = new ArrayList<>();
                    currentFeatureList.add(f);
                    currentMin = start;
                    currentMax = end;

                } else if (currentFeatureList.size() > batchSize) {

                    // New interval
                    IntervalTree<List<Feature>> tree = featureMap.get(lastChr);
                    tree.insert(new Interval(currentMin, currentMax, currentFeatureList));

                    currentFeatureList = new ArrayList<>();
                    currentFeatureList.add(f);
                    currentMin = start;
                    currentMax = end;


                } else {

                    // Update interval
                    currentMin = Math.min(currentMin, start);
                    currentMax = Math.max(currentMax, end);
                    currentFeatureList.add(f);

                }
            }
        }
    }

    private static Comparator<Feature> getPositionComparator() {
        Comparator<Feature> comp = new Comparator<Feature>() {
            private Comparator<String> nameComparator = ChromosomeNameComparator.get();
            public int compare(Feature o1, Feature o2) {
                int nameComp = nameComparator.compare(o1.getChr(), o2.getChr());
                if (nameComp == 0) {
                    return o1.getStart() - o2.getStart();
                } else {
                    return nameComp;
                }
            }
        };
        return comp;
    }
}
