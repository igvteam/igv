/*
 * FeatureUtils.java
 *
 * Useful utilities for working with Features
 */
package org.igv.feature;

import com.google.common.base.Predicate;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;

import java.util.*;

/**
 * @author jrobinso
 */
public class FeatureUtils {


    public static Predicate<Feature> getOverlapPredicate(final String chr, final int start, final int end) {
        return object -> chr.equals(object.getChr()) && object.getStart() <= end && object.getEnd() > start;
    }


    /**
     * Sort the feature list by ascending start value
     */
    public static void sortFeatureList(List<? extends Locatable> features) {
        Collections.sort(features, FEATURE_START_COMPARATOR);
    }

    /**
     * Return a feature from the supplied list whose extent, expanded by "buffer", contains the given position.
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
     * Return the first feature whose start is > position
     *
     * @param position
     * @param features
     * @return
     */
    public static Feature getFeatureStartsAfter(double position, List<? extends Feature> features) {
        if (features.size() == 0 ||
                features.get(features.size() - 1).getStart() <= position) {
            return null;
        }
        int idxBefore = getIndexBefore(position, features);
        if (idxBefore >= features.size() - 1) {
            return null;
        } else {
            for (Feature f : features) {
                if (f.getStart() > position) return f;
            }
        }
        return null;
    }

    /**
     * Return the first feature whose center is > the given position.  If no features satisfy the criteria
     * return null;
     *
     * @param position
     * @param features
     * @return
     */
    public static Feature getFeatureCenteredAfter(double position, List<? extends Feature> features) {

        if (features.size() == 0 ||
                center(features.get(features.size() - 1)) <= position) {
            return null;
        }

        int idx = getIndexCenterAfter(position, features);
        if (idx < 0 || idx > features.size() - 1) {
            return null;
        } else {
            return features.get(idx);
        }
    }

    /**
     * Return the first feature whose center is < the given position.  If no features satisfy the criteria
     * return null;
     *
     * @param position
     * @param features
     * @return
     */
    public static Feature getFeatureCenteredBefore(double position, List<? extends Feature> features) {

        if (features.size() == 0) {
            return null;
        }
        if (center(features.get(0)) >= position) {
            return null;
        }


        int idx = getIndexCenterAfter(position, features);
        if (idx < 0) {
            return null;
        } else {
            idx--;
            while (idx >= 0) {
                if (center(features.get(idx)) < position) {
                    return features.get(idx);
                }
                idx--;
            }
        }
        return null;
    }

    /**
     * Get the feature which is closest to the given position.
     * If multiple features are equidistant it breaks ties by the one that has the closest
     * starting position, followed by the shortest feature.
     */
    public static Feature getFeatureClosest(double position, List<? extends Feature> features) {
        final Comparator<Feature> closestComparator = Comparator.comparing((Feature f) -> getDistance((int) position, f)) //overlaps or not
                .thenComparing(f -> Math.abs(position - f.getStart()))
                .thenComparing(Feature::getLengthOnReference);
        return Collections.min(features, closestComparator);
    }

    private static int getDistance(final int position, final Feature f) {
        final int end = f.getEnd();
        final int start = f.getStart();
        if (position >= start && position < end) { //within the feature
            return 0;
        } else if (position >= end) {
            return position - end + 1; //after the feature
        } else {
            return start - position; //before the feature
        }
    }

    /**
     * Return the index to the last feature in the list with a start < the given position.  It is assumed
     * the list is sorted by start position.
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
            return -1;
        }

        KeyClass key = new KeyClass(position);
        int idx = Collections.binarySearch(features, key, FEATURE_START_COMPARATOR);
        if (idx < 0) idx = -1 * idx;
        idx = Math.min(features.size() - 1, idx);
        while (idx > 0) {
            if (features.get(idx).getStart() < position) {
                break;
            } else {
                idx--;
            }
        }
        return idx;
    }

    /**
     * Return the index to the first feature in the list with a center > the given position.  It is assumed
     * the list is sorted by center position.
     *
     * @param position
     * @param features
     * @return
     */
    private static int getIndexCenterAfter(double position, List<? extends Feature> features) {

        Feature first = features.get(0);
        Feature last = features.get(features.size() - 1);
        if (center(first) > position) {
            return 0;
        }

        KeyClass key = new KeyClass(position);
        int idx = Collections.binarySearch(features, key, FEATURE_CENTER_COMPARATOR);
        if (idx < 0) idx = -1 * idx - 1;
        return idx;
    }

    /**
     * Return the index to the first feature in the list with a center < the given position.  It is assumed
     * the list is sorted by center position.  If no features satisfies the criteria return -1
     *
     * @param position
     * @param features
     * @return
     */
    public static int getIndexCenterBefore(double position, List<? extends Feature> features) {

        if (features == null || features.size() == 0) {
            return -1;
        }

        Feature first = features.get(0);
        Feature last = features.get(features.size() - 1);
        if (center(first) < position) {
            return 0;
        }
        if (center(last) >= position) {
            return features.size();
        }

        KeyClass key = new KeyClass(position);
        int idx = Collections.binarySearch(features, key, FEATURE_CENTER_COMPARATOR);
        if (idx < 0) idx = -1 * idx;
        idx = Math.min(features.size() - 1, idx);
        return idx;
    }

    private static double center(Feature f) {
        return (f.getStart() + f.getEnd()) / 2.0;
    }

    /**
     * Return the index to the last feature in the list with a start < the given position
     *
     * @param position
     * @param features
     * @return
     */
    public static int getIndexBeforeOld(double position, List<? extends Feature> features) {

        if (features == null || features.size() == 0) {
            return -1;
        }
        if (features.get(features.size() - 1).getStart() <= position) {
            return features.size() - 1;
        }
        if (features.get(0).getStart() >= position) {
            return -1;
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
     * Return all features from the supplied list who's extent, expanded to "minWidth" if needed,  contains the given position
     *
     * @param position
     * @param flanking - flanking distance in BP for hit test,  a hit occurs over the feature region +/- flanking
     * @param features
     * @return
     */
    public static List<Feature> getAllFeaturesContaining(double position,
                                                         double flanking,
                                                         List<? extends htsjdk.tribble.Feature> features) {

        List<Feature> returnList = null;

        for (Feature feature : features) {
            if (feature.getStart() > position + flanking) {
                break;
            }
            if (feature.getEnd() < position - flanking) {
                continue;
            }
            if (returnList == null) returnList = new ArrayList();
            returnList.add(feature);
        }

        // Sort features by distance from position (features closest to position listed first)
        if (returnList != null) {
            returnList.sort((o1, o2) -> {
                double dist1 = Math.abs((o1.getStart() + (o1.getEnd() - o1.getStart()) / 2) - position);
                double dist2 = Math.abs((o2.getStart() + (o2.getEnd() - o2.getStart()) / 2) - position);
                return dist1 == dist2 ? 0 : dist1 > dist2 ? 1 : -1;
            });
        }
        return returnList;
    }

    public static final Comparator<Locatable> FEATURE_START_COMPARATOR = (o1, o2) -> o1.getStart() - o2.getStart();
    public static final Comparator<Feature> FEATURE_END_COMPARATOR = (o1, o2) -> o1.getEnd() - o2.getEnd();
    public static final Comparator<Feature> FEATURE_CENTER_COMPARATOR = (o1, o2) -> o1.getStart() - o2.getStart() + o1.getEnd() - o2.getEnd();

    /**
     * Compute reading frames
     *
     * @param gene
     */
    public static void computeReadingFrames(IGVFeature gene) {

        List<Exon> exons = gene.getExons();
        if (exons.size() == 0) {
            return;
        }

        int startIndex = (gene.getStrand() == Strand.POSITIVE) ? 0 : exons.size() - 1;
        int endIndex = (gene.getStrand() == Strand.POSITIVE) ? exons.size() : -1;
        int increment = (gene.getStrand() == Strand.POSITIVE) ? 1 : -1;
        int cds = 0;
        int exonNumber = 1;
        for (int i = startIndex; i != endIndex; i += increment) {

            Exon exon = exons.get(i);
            exon.setNumber(exonNumber++);

            if (exon.getCodingLength() > 0 || cds > 0) {  // Skip until we find the coding start
                int modCds = cds % 3;
                int frame = modCds;  //(modCds == 0) ? 0 : 3 - modCds;
                exon.setReadingFrame(frame);
                cds += exon.getCodingLength();
            }
        }
    }
}

class KeyClass implements Feature {

    double position;

    public KeyClass(double position) {
        this.position = position;
    }

    @Override
    public String getContig() {
        return null;
    }

    @Override
    public int getStart() {
        return (int) (position);
    }

    @Override
    public int getEnd() {
        return (int) position + 1;
    }
}
