/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

/*
 * FeatureUtils.java
 *
 * Useful utilities for working with Features
 */
package org.broad.igv.feature;

import com.google.common.base.Predicate;
import org.broad.tribble.Feature;

import java.util.*;

/**
 * @author jrobinso
 */
public class FeatureUtils {


    public static Predicate<Feature> getOverlapPredicate(final String chr, final int start, final int end) {
        Predicate<Feature> overlapPredicate = new Predicate<Feature>() {
            @Override
            public boolean apply(Feature object) {
                return chr.equals(object.getChr()) && object.getStart() <= end && object.getEnd() >= start;
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
    public static List<List<IGVFeature>> segreateFeatures(List<IGVFeature> features, double scale) {

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
    public static Feature getFeatureAt(double position, int buffer, List<? extends Feature> features) {

        int startIdx = 0;
        int endIdx = features.size();

        while (startIdx != endIdx) {
            int idx = (startIdx + endIdx) / 2;

            org.broad.tribble.Feature feature = features.get(idx);

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
            org.broad.tribble.Feature f = features.get(index);
            if (f.getStart() < position) {
                return f;
            }
            index--;
        }
        return null;

    }

    public static Feature getFeatureClosest(double position, List<? extends org.broad.tribble.Feature> features) {
        // look for exact match at position:
        org.broad.tribble.Feature f0 = getFeatureAt(position, features);
        if (f0 != null) {
            return f0;
        }
        // otherwise look for features on either side and return the closest:
        org.broad.tribble.Feature f1 = getFeatureBefore(position, features);
        org.broad.tribble.Feature f2 = getFeatureAfter(position, features);

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
                                                 List<? extends org.broad.tribble.Feature> features) {

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
