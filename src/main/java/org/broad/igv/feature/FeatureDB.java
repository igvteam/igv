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

package org.broad.igv.feature;

//~--- non-JDK imports --------------------------------------------------------

import com.jidesoft.utils.SortedList;
import htsjdk.tribble.Feature;
import htsjdk.tribble.NamedFeature;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;

import java.util.*;

/**
 * Used to support the "search" box.
 *
 * @author jrobinso
 */
public class FeatureDB {

    private static Logger log = LogManager.getLogger(FeatureDB.class);

    private Map<String, List<NamedFeature>> featureMap = Collections.synchronizedSortedMap(new TreeMap<>());
    private final int MAX_DUPLICATE_COUNT = 20;

    public void addFeature(NamedFeature feature) {

        final String name = feature.getName();
        if (name != null && name.length() > 0 && !name.equals(".")) {
            put(name, feature);
        }
        if (feature instanceof BasicFeature) {
            final BasicFeature igvFeature = (BasicFeature) feature;
            final String id = igvFeature.getIdentifier();
            if (id != null && id.length() > 0) {
                put(id, feature);
            }

            addByAttributes(igvFeature);

            List<Exon> exons = igvFeature.getExons();
            if (exons != null) {
                for (Exon exon : exons) {
                    addByAttributes(exon);
                }
            }
        }
    }

    private void addByAttributes(IGVFeature igvFeature) {
        List<String> attributeKeys = igvFeature.getAttributeKeys();
        for (String key : attributeKeys) {
            String value = igvFeature.getAttribute(key);
            if (value.length() < 50) {
                put(value, igvFeature);
            }
        }
    }

    /**
     * Add feature to the list of features associated with this name.
     * Performs no data integrity checks
     *
     * @param name
     * @param feature
     * @return true if successfully added, false if not
     */
    void put(String name, NamedFeature feature) {

        String key = name.toUpperCase();

        // Use computeIfAbsent for safer and more efficient concurrent access
        List<NamedFeature> currentList = featureMap.computeIfAbsent(key, k ->
                new SortedList<>(new ArrayList<>(), FeatureComparator.get(true)));

        // The list itself is not thread-safe, so synchronize on it during modification.
        synchronized (currentList) {
            if (currentList.size() < MAX_DUPLICATE_COUNT) {
                currentList.add(feature);
            }
        }
    }


    public void addFeature(String name, IGVNamedFeature feature) {
        put(name.toUpperCase(), feature);
    }

    public void addFeatures(List<htsjdk.tribble.Feature> features) {
        for (htsjdk.tribble.Feature feature : features) {
            if (feature instanceof IGVFeature)
                addFeature((IGVFeature) feature);
        }
    }


    public void clearFeatures() {
        featureMap.clear();
    }

    int size() {
        return featureMap.size();
    }

    /**
     * Return a feature with the given name.
     */
    public NamedFeature getFeature(String name) {
        String nm = name.trim().toUpperCase();
        List<NamedFeature> features = featureMap.get(nm);

        if (features != null) {
            return features.get(0);
        } else {
            return null;
        }
    }

    /**
     * Return all features with an exact match to the given name.
     */
    public List<NamedFeature> getFeaturesMatching(String name) {
        String nm = name.trim().toUpperCase();
        return featureMap.get(nm);
    }

    /**
     * Get all features which match nm. Not necessarily an exact match. Current implementation will match anything
     * for which name is at the beginning.
     * Note: This method is synchronized on featureMap
     * @param name : Search string. Features which begin with this string will be found.
     * @return
     */
    Map<String, List<NamedFeature>> getFeaturesMap(String name) {
        String nm = name.trim().toUpperCase();
        SortedMap<String, List<NamedFeature>> treeMap = (SortedMap) featureMap;
        //Search is inclusive to first argument, exclusive to second
        return treeMap.subMap(nm, nm + Character.MAX_VALUE);
    }

    /**
     * Get a list of features which start with the provided name.  In the case of multiple features with the same
     * name only one will be returned.  This is used for the text hints in the search box.  Value questionable.
     *
     * Note that matches can be inexact
     *
     * @param name
     * @param limit
     * @return
     */
    public List<NamedFeature> getFeaturesStartingWith(String name, int limit) {

        //Note: We are iterating over submap, this needs
        //to be synchronized over the main map.
        synchronized (featureMap) {
            Map<String, List<NamedFeature>> resultMap = getFeaturesMap(name);
            Set<String> names = resultMap.keySet();
            Iterator<String> nameIter = names.iterator();
            ArrayList<NamedFeature> features = new ArrayList<>((Math.min(limit, names.size())));
            int ii = 0;
            while (nameIter.hasNext() && ii < limit) {
                List<NamedFeature> subFeats = resultMap.get(nameIter.next());
                features.add(subFeats.get(0));
                ii++;
            }
            return features;
        }
    }

    /**
     * Doubleton class. Can sort forward or descending, at most 2 instances.
     */
    private static class FeatureComparator implements Comparator<Feature> {
        private boolean descending;
        private static FeatureComparator ascending_instance;
        private static FeatureComparator descending_instance;

        public static FeatureComparator get(boolean descending) {
            FeatureComparator instance;
            if (descending) {
                if (ascending_instance == null) {
                    ascending_instance = new FeatureComparator(descending);
                }
                instance = ascending_instance;
            } else {
                if (descending_instance == null) {
                    descending_instance = new FeatureComparator(descending);
                }
                instance = descending_instance;
            }

            return instance;
        }

        private FeatureComparator(boolean reverse) {
            this.descending = reverse;
        }

        public int compare(Feature feat1, Feature feat2) {

            // Prefer the shortest chromosome name.  Longer names are most likely "weird"
            // e.g.  chr1_gl000191_random
            int nameLen1 = feat1.getChr().length();
            int nameLen2 = feat2.getChr().length();
            if (nameLen1 != nameLen2) {
                return nameLen1 - nameLen2;
            }


            int len1 = (feat1.getEnd() - feat1.getStart());
            int len2 = (feat2.getEnd() - feat2.getStart());
            int toRet;
            if (!this.descending) {
                toRet = len1 - len2;
            } else {
                toRet = len2 - len1;
            }

            return toRet;

        }
    }

}
