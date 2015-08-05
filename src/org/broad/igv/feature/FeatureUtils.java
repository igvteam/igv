/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/*
 * FeatureUtils.java
 *
 * Useful utilities for working with Features
 */
package org.broad.igv.feature;

import com.google.common.base.Predicate;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import htsjdk.tribble.Feature;

import java.util.*;

/**
 * @author jrobinso
 */
public class FeatureUtils {


    public static Predicate<Feature> getOverlapPredicate(final String chr, final int start, final int end) {
        Predicate<Feature> overlapPredicate = new Predicate<Feature>() {
            @Override
            public boolean apply(Feature object) {
                return chr.equals(object.getChr()) && object.getStart() <= end && object.getEnd() > start;
            }
        };
        return overlapPredicate;
    }

    public static Map<String, List<IGVFeature>> divideByChromosome(List<IGVFeature> features) {
        Map<String, List<IGVFeature>> featureMap = new LinkedHashMap();
        for (IGVFeature f : features) {
            List<IGVFeature> flist = featureMap.get(f.getChr());
            if (flist == null) {
                flist = new ArrayList();
                featureMap.put(f.getChr(), flist);
            }
            flist.add(f);
        }
        return featureMap;
    }

    /**
     * Segregate a list of possibly overlapping features into a list of
     * non-overlapping lists of features.
     */
    public static List<List<IGVFeature>> segregateFeatures(List<IGVFeature> features, double scale) {

        // Create a list to hold the lists of non-overlapping features
        List<List<IGVFeature>> segmentedLists = new ArrayList();

        // Make a working copy of the original list.
        List<IGVFeature> workingList = new LinkedList(features);
        sortFeatureList(workingList);

        // Loop until all features have been allocated to non-overlapping lists
        while (workingList.size() > 0) {

            List<IGVFeature> nonOverlappingFeatures = new LinkedList();
            List<IGVFeature> overlappingFeatures = new LinkedList();

            // Prime the loop with the first feature, it can't overlap itself
            IGVFeature f1 = workingList.remove(0);
            nonOverlappingFeatures.add(f1);
            while (workingList.size() > 0) {
                IGVFeature f2 = workingList.remove(0);
                int scaledStart = (int) (f2.getStart() / scale);
                int scaledEnd = (int) (f1.getEnd() / scale);
                if (scaledStart > scaledEnd) {
                    nonOverlappingFeatures.add(f2);
                    f1 = f2;
                } else {
                    overlappingFeatures.add(f2);
                }
            }

            // Add the list of non-overlapping features and start again with whats left
            segmentedLists.add(nonOverlappingFeatures);
            workingList = overlappingFeatures;
        }
        return segmentedLists;
    }

    /**
     * Sort the feature list by ascending start value
     */
    public static void sortFeatureList(List<? extends Feature> features) {
        Collections.sort(features, FEATURE_START_COMPARATOR);
    }


    /**
     * Null safe version of {@linkplain #combineSortedFeatureListsNoDups(java.util.Iterator, java.util.Iterator, int, int)}
     * If BOTH self and other are null, returns null. If only one is null,
     * returns the other
     *
     * @param self
     * @param other
     * @param start
     * @param end
     * @return
     */
    public static <T extends Feature> List<T> combineSortedFeatureListsNoDups(List<T> self, List<T> other, int start, int end) {
        if (self == null && other == null) {
            return null;
        } else if (self == null) {
            return other;
        } else if (other == null) {
            return self;
        }

        return combineSortedFeatureListsNoDups(self.iterator(), other.iterator(), start, end);
    }

    /**
     * Features are sorted by start position. The interval being merged
     * will have some features on the left or right that the current
     * interval does not have. Both are sorted by start position.
     * So we first add at the beginning, and then the end,
     * only those features which don't overlap the original interval.
     *
     * @param selfIter  iterator of features belonging to this interval
     * @param otherIter iterator of features belonging to some other interval
     * @param start     the beginning of the interval from which selfIter was derived
     * @param end       the end of the interval from which selfIter was derived
     * @return Combined sorted list.
     * @throws ClassCastException If the elements of an iterator cannot be cast
     *                            to a Feature.
     */
    public static <T extends Feature> List<T> combineSortedFeatureListsNoDups(Iterator<T> selfIter, Iterator<T> otherIter, int start, int end) {
        List<T> allFeatures = new ArrayList<T>();
        T otherFeat = null;

        while (otherIter.hasNext()) {
            otherFeat = otherIter.next();
            if (otherFeat.getEnd() > start) break;
            allFeatures.add(otherFeat);
        }

        while (selfIter.hasNext()) {
            allFeatures.add(selfIter.next());
        }

        while (otherIter.hasNext()) {
            if (otherFeat.getStart() >= end) {
                allFeatures.add(otherFeat);
            }
            otherFeat = otherIter.next();
        }

        if (otherFeat != null && otherFeat.getStart() >= end) {
            allFeatures.add(otherFeat);
        }

        return allFeatures;
    }

    /**
     * Return a feature from the supplied list at the given position.
     *
     * @param position 0-based genomic position to which to search for feature
     * @param buffer   search region. The first feature which contains the start position, (expanded by buffer, inclusive)
     *                 will be accepted.
     * @param features
     * @return
     */
    public static <T extends Feature> T getFeatureAt(double position, int buffer, List<? extends T> features) {

        int startIdx = 0;
        int endIdx = features.size();

        while (startIdx != endIdx) {
            int idx = (startIdx + endIdx) / 2;

            T feature = features.get(idx);

            int effectiveStart = feature.getStart();
            int effectiveEnd = feature.getEnd();

            if (position >= effectiveStart - buffer) {
                if (position <= effectiveEnd + buffer) {
                    return features.get(idx);
                } else {
                    if (idx == startIdx) {
                        return null;
                    } else {
                        startIdx = idx;
                    }
                }
            } else {
                endIdx = idx;
            }
        }

        return null;
    }

    /**
     * Get the index of the feature just to the right of the given position.
     * If there is no feature to the right return -1;
     *
     * @param position
     * @param features
     * @return
     */
    public static Feature getFeatureAfter(double position, List<? extends Feature> features) {

        if (features.size() == 0 ||
                features.get(features.size() - 1).getStart() <= position) {
            return null;
        }

        int startIdx = 0;
        int endIdx = features.size();

        // Narrow the list to ~ 10
        while (startIdx != endIdx) {
            int idx = (startIdx + endIdx) / 2;
            double distance = features.get(idx).getStart() - position;
            if (distance <= 0) {
                startIdx = idx;
            } else {
                endIdx = idx;
            }
            if (endIdx - startIdx < 10) {
                break;
            }
        }

        // Now find feature
        for (int idx = startIdx; idx < features.size(); idx++) {
            if (features.get(idx).getStart() > position) {
                return features.get(idx);
            }
        }

        return null;

    }

    public static Feature getFeatureBefore(double position, List<? extends Feature> features) {

        int index = getIndexBefore(position, features);
        while (index >= 0) {
            htsjdk.tribble.Feature f = features.get(index);
            if (f.getStart() < position) {
                return f;
            }
            index--;
        }
        return null;

    }

    public static Feature getFeatureClosest(double position, List<? extends htsjdk.tribble.Feature> features) {
        // look for exact match at position:
        htsjdk.tribble.Feature f0 = getFeatureAt(position, features);
        if (f0 != null) {
            return f0;
        }
        // otherwise look for features on either side and return the closest:
        htsjdk.tribble.Feature f1 = getFeatureBefore(position, features);
        htsjdk.tribble.Feature f2 = getFeatureAfter(position, features);

        double d1 = f1 == null ? Double.MAX_VALUE : Math.abs(position - f1.getEnd());
        double d2 = f2 == null ? Double.MAX_VALUE : Math.abs(f2.getStart() - position);

        return (d1 < d2 ? f1 : f2);

    }

    /**
     * Return a feature that encompasses the supplied position.
     *
     * @param position Query position.
     * @param features List of features.
     * @return The feature whose start overlaps with position, or null.
     */
    private static Feature getFeatureAt(double position, List<? extends Feature> features) {
        int strt = (int) position;
        Feature key = new BasicFeature("", strt, strt + 1);

        int r = Collections.binarySearch(features, key, FEATURE_START_COMPARATOR);

        if (r >= 0) {
            return features.get(r);
        } else {
            return null;
        }
    }


    /**
     * Return the index to the last feature in the list with a start < the given position
     *
     * @param position
     * @param features
     * @return
     */
    public static int getIndexBefore(double position, List<? extends Feature> features) {

        if (features == null || features.size() == 0) {
            return -1;
        }
        if (features.get(features.size() - 1).getStart() <= position) {
            return features.size() - 1;
        }
        if (features.get(0).getStart() >= position) {
            return 0;
        }

        int startIdx = 0;
        int endIdx = features.size() - 1;

        while (startIdx != endIdx) {
            int idx = (startIdx + endIdx) / 2;
            double distance = features.get(idx).getStart() - position;
            if (distance <= 0) {
                startIdx = idx;
            } else {
                endIdx = idx;
            }
            if (endIdx - startIdx < 10) {
                break;
            }
        }

        if (features.get(endIdx).getStart() >= position) {
            for (int idx = endIdx; idx >= 0; idx--) {
                if (features.get(idx).getStart() < position) {
                    return idx;
                }
            }
        } else {
            for (int idx = endIdx + 1; idx < features.size(); idx++) {
                if (features.get(idx).getStart() >= position) {
                    return idx - 1;
                }

            }
        }
        return -1;
    }

    /**
     * Return a feature from the supplied list at the given position.
     *
     * @param position
     * @param maxLength
     * @param features
     * @return
     */
    public static List<Feature> getAllFeaturesAt(double position,
                                                 double maxLength,
                                                 double minWidth,
                                                 List<? extends htsjdk.tribble.Feature> features) {

        List<Feature> returnList = null;

        double adjustedPosition = Math.max(0, position - maxLength);
        int startIdx = Math.max(0, getIndexBefore(adjustedPosition, features));
        for (int idx = startIdx; idx < features.size(); idx++) {
            Feature feature = features.get(idx);
            int start = feature.getStart() - (int) (minWidth / 2);

            if (start > position) {
                break;
            }

            int end = feature.getEnd() + (int) (minWidth / 2);

            if (position >= start && position <= end) {
                if (returnList == null) returnList = new ArrayList();
                returnList.add(feature);
            }
        }

        return returnList;
    }


    /**
     * Export a gene in "edx" format
     * LYZ
     * hg19 12	+	69742133	69748013	69742188	69746999
     * 4
     * 69742133	69742324	0
     * 69743887	69744052	1
     * 69745999	69746078	1
     * 69746932	69748013	2
     * 69742083	69742374	AAGGGTGGAGCC
     *
     * @param feature
     */
    public static void exportEDX(IGVFeature feature) {

        Genome genome = GenomeManager.getInstance().getCurrentGenome();

        System.out.println(feature.getName());
        System.out.print(genome.getId() + "\t");
        String str = feature.getStrand() == Strand.POSITIVE ? "+" : "-";
        System.out.print(feature.getChr() + "\t" + str + "\t" + feature.getStart() + "\t" + feature.getEnd() + "\t");

        List<Exon> exons = feature.getExons();

        // Find coding start & end
        int cdStart = 0;
        int cdEnd = 0;
        for (Exon ex : exons) {
            if (ex.getCdStart() > ex.getStart()) {
                cdStart = ex.getCdStart();
                break;
            }
        }
        for (int i = 0; i < exons.size(); i++) {
            Exon ex = exons.get(i);
            if (ex.getCdEnd() < ex.getEnd()) {
                cdEnd = ex.getCdEnd();
                break;
            }
        }
        System.out.println(cdStart + "\t" + cdEnd);
        System.out.println(exons.size());

        for(Exon ex : exons) {
            int rs = ex.getReadingFrame();
            int fs = rs == 0 ? 0 : 3 - rs;
            System.out.println(ex.getStart() + "\t" + ex.getEnd() + "\t" + fs);
        }


        // Sequence
        int buffer = 10;
        for(int i=0; i<exons.size(); i++) {
            Exon ex = exons.get(i);
            // Buffer -- 50 bp by default but can be smaller

            int seqStart = ex.getStart() - buffer;

            int nextExonStart = (i < exons.size() - 1 ? exons.get(i+1).getStart() : Integer.MAX_VALUE);
            buffer = Math.min(50, (nextExonStart - ex.getEnd()) / 2);
            int seqEnd = ex.getEnd() + buffer;

            byte [] sequence = genome.getSequence(feature.getChr(), seqStart, seqEnd);
            String seqString = new String(sequence);
            System.out.println(seqStart + "\t" + seqEnd + "\t" + seqString);
        }

    }


    private static final Comparator<Feature> FEATURE_CONTAINS_COMPARATOR = new Comparator<Feature>() {
        public int compare(Feature o1, Feature o2) {
            int genomeStart2 = o2.getStart();
            int genomeStart1 = o1.getEnd();
            if (genomeStart2 >= genomeStart1 && o2.getEnd() <= o1.getEnd()) {
                return 0;
            } else {
                return genomeStart1 - genomeStart2;
            }
        }
    };

    public static final Comparator<Feature> FEATURE_START_COMPARATOR = new Comparator<Feature>() {
        public int compare(Feature o1, Feature o2) {
            return o1.getStart() - o2.getStart();
        }
    };
}
