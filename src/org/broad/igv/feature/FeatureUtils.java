/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

/*
 * FeatureUtils.java
 *
 * Useful utilities for working with Features
 */
package org.broad.igv.feature;

import org.broad.tribble.Feature;

import java.util.*;

/**
 * @author jrobinso
 */
public class FeatureUtils {

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

        Collections.sort(features, new Comparator() {

            public int compare(Object o1, Object o2) {
                org.broad.tribble.Feature f1 = (org.broad.tribble.Feature) o1;
                org.broad.tribble.Feature f2 = (org.broad.tribble.Feature) o2;
                return (f1.getStart() - f2.getStart());
            }
        });
    }

    /**
     * Return a feature from the supplied list at the given position.
     *
     * @param position
     * @param buffer
     * @param features
     * @return
     */
    public static Feature getFeatureAt(double position, int buffer, List<? extends Feature> features) {

        position--;

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

        double d1 = f1 == null ? Double.MAX_VALUE : Math.abs(f1.getEnd() - position);
        double d2 = f2 == null ? Double.MAX_VALUE : Math.abs(f2.getStart() - position);

        return (d1 < d2 ? f1 : f2);

    }

    /**
     * Return a feature whose position is exactly 'position'.
     * @param position Query position.
     * @param features List of features.
     * @return The feature whose start overlaps with position, or null.
     */
    private static Feature getFeatureAt(double position, List<? extends Feature> features) {
        Feature key = new BasicFeature("", (int) position, (int) position);
        int r = Collections.binarySearch(features, key, new Comparator<Object>() {
            public int compare(Object o1, Object o2) {
                Feature f1 = (Feature) o1;
                Feature f2 = (Feature) o2;
                return f1.getEnd() - f2.getEnd();
            }
        });

        if (r > 0) {
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
}
