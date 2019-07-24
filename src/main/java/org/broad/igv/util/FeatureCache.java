package org.broad.igv.util;

import htsjdk.samtools.util.Locatable;
import org.broad.igv.feature.genome.ChromosomeNameComparator;

import java.util.*;

public class FeatureCache<T extends Locatable> {

    int DEFAULT_BATCH_SIZE = 100;

    int batchSize;
    Map<String, IntervalTree<List<T>>> featureMap;

    public FeatureCache(List<T> features) {
        init(features, DEFAULT_BATCH_SIZE);
    }

    public FeatureCache(List<T> features, int batchSize) {
        init(features, batchSize);
    }

    public List<T> getFeatures(String chr, int start, int end) {
        List<T> features = new ArrayList<>();
        IntervalTree<List<T>> tree = featureMap.get(chr);
        if(tree != null) {
            List<Interval<List<T>>> intervals = tree.findOverlapping(start, end);
            for (Interval<List<T>> interval : intervals) {
                List<T> intervalFeatures = interval.getValue();
                for (T f : intervalFeatures) {
                    if (f.getEnd() >= start && f.getStart() <= end) {
                        features.add(f);
                    }
                }
            }
        }
        return features;
    }

    private void init(List<T> features, int batchSize) {
        String lastChr = null;

        List<T> currentFeatureList = new ArrayList<>();
        int currentMin = Integer.MAX_VALUE;
        int currentMax = 0;

        featureMap = new HashMap<>();

        features.sort(getPositionComparator());

        for (T f : features) {

            final String chr = f.getContig();
            final int start = f.getStart();
            final int end = f.getEnd();

            if (lastChr == null) {
                currentMin = start;
                currentMax = end;
                currentFeatureList.add(f);
                IntervalTree<List<T>> tree = new IntervalTree<>();
                featureMap.put(chr, tree);
                lastChr = chr;
            } else {

                if (!chr.equals(lastChr)) {

                    // New tree
                    IntervalTree<List<T>> tree = featureMap.get(lastChr);
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
                    IntervalTree<List<T>> tree = featureMap.get(lastChr);
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

        IntervalTree<List<T>> tree = featureMap.get(lastChr);
        tree.insert(new Interval(currentMin, currentMax, currentFeatureList));


    }

    private  Comparator<T> getPositionComparator() {
        Comparator<T> comp = new Comparator<T>() {
            private Comparator<String> nameComparator = ChromosomeNameComparator.get();
            public int compare(T o1, T o2) {
                int nameComp = nameComparator.compare(o1.getContig(), o2.getContig());
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
