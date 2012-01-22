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
package org.broad.igv.feature;

//~--- non-JDK imports --------------------------------------------------------

import com.jidesoft.utils.SortedList;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ui.IGV;
import org.broad.tribble.Feature;

import java.util.*;

/**
 * This is a placeholder class for a true "feature database" wrapper.  Its purpose
 * is to return a feature given a name.  Used to support the "search" box.
 *
 * @author jrobinso
 */
public class FeatureDB {

    private static Logger log = Logger.getLogger(FeatureDB.class);
    /**
     * Map for all features other than genes.
     */
    //private static Map<String, NamedFeature> featureMap = new HashMap(10000);
    private static Map<String, List<NamedFeature>> featureMap = Collections.synchronizedSortedMap(new TreeMap<String, List<NamedFeature>>());
    private static Genome GENOME = null;
    private static final int MAX_DUPLICATE_COUNT = 20;

    public static void addFeature(NamedFeature feature) {

        final String name = feature.getName();
        if (name != null && name.length() > 0 && !name.equals(".")) {
            put(name, feature);
        }
        if (feature instanceof IGVFeature) {
            final IGVFeature igvFeature = (IGVFeature) feature;
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

    private static void addByAttributes(IGVFeature igvFeature) {
        Map<String, String> attributes = igvFeature.getAttributes();
        if (attributes != null) {
            for (String value : attributes.values()) {
                if (value.length() < 20) {
                    put(value, igvFeature);
                }
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
    public static boolean put(String name, NamedFeature feature) {
        String key = name.toUpperCase();
        Genome currentGenome = GENOME;
        if (!Globals.isHeadless()) {
            currentGenome = IGV.getInstance().getGenomeManager().getCurrentGenome();
        }
        if (currentGenome != null && currentGenome.getChromosome(feature.getChr()) == null) {
            return false;
        }

        synchronized (featureMap) {
            List<NamedFeature> currentList = featureMap.get(key);
            if (currentList == null) {
                List<NamedFeature> newList = new SortedList<NamedFeature>(
                        new ArrayList<NamedFeature>(), FeatureComparator.get(true));
                boolean added = newList.add(feature);
                if (added) {
                    featureMap.put(key, newList);
                }
                return added;
            } else {
                // Don't let list grow without bounds
                if (currentList.size() > MAX_DUPLICATE_COUNT) {
                    return false;
                }

                //Don't want to add exact duplicates. JTR -- why not?
//                for (NamedFeature f : currentList) {
//                    if ((f.getStart() == feature.getStart()) && (f.getEnd() == feature.getEnd())) {
//                        return false;
//                    }
//                }
                return currentList.add(feature);
            }
        }
    }

    static void setGenome(Genome genome) {
        GENOME = genome;
    }


    public static void addFeature(String name, NamedFeature feature) {
        put(name.toUpperCase(), feature);
    }


    private FeatureDB() {
    }


    public static void addFeatures(List<org.broad.tribble.Feature> features) {
        for (org.broad.tribble.Feature feature : features) {
            if (feature instanceof IGVFeature)
                addFeature((IGVFeature) feature);
        }
    }


    public static void clearFeatures() {
        featureMap.clear();
    }

    static int size() {
        return featureMap.size();
    }

    /**
     * Return the feature, if any, with the given name.  Genes are given
     * precedence.
     */
    public static NamedFeature getFeature(String name) {
        String nm = name.trim().toUpperCase();
        List<NamedFeature> features = featureMap.get(nm);

        if (features != null) {
            return features.get(0);
        } else {
            return null;
        }
    }

    /**
     * Get all features which match nm. Not necessarily
     * an exact match. Current implementation will match anything
     * for which name is at the beginning, including but not limited to
     * exact matches.
     * <p/>
     * NOTE: "It is imperative that the user manually synchronize
     * on [this sorted map] when iterating over any of its
     * collection views, or the collections views of any of its
     * subMap, headMap or tailMap views". See
     * <a href="http://docs.oracle.com/javase/6/docs/api/java/util/Collections.html#synchronizedSortedMap%28java.util.SortedMap%29> here</a>
     *
     * @param name : Search string. Features which begin with this
     *             string will be found.
     * @return
     */
    static Map<String, List<NamedFeature>> getFeaturesMap(String name) {
        String nm = name.trim().toUpperCase();
        SortedMap<String, List<NamedFeature>> treeMap = (SortedMap) featureMap;
        //Search is inclusive to first argument, exclusive to second
        return treeMap.subMap(nm, nm + Character.MAX_VALUE);
    }

    public static List<NamedFeature> getFeaturesList(String name, int limit) {
        return getFeaturesList(name, limit, true);
    }

    public static List<NamedFeature> getFeaturesList(String name, int limit, boolean longestOnly) {

        //Note: We are iterating over submap, this needs
        //to be synchronized over the main map.
        synchronized (featureMap) {
            Map<String, List<NamedFeature>> resultMap = getFeaturesMap(name);
            Set<String> names = resultMap.keySet();
            Iterator<String> nameIter = names.iterator();
            ArrayList<NamedFeature> features = new ArrayList<NamedFeature>((Math.min(limit, names.size())));
            int ii = 0;
            while (nameIter.hasNext() && ii < limit) {
                List<NamedFeature> subFeats = resultMap.get(nameIter.next());
                if (longestOnly) {
                    features.add(subFeats.get(0));
                } else {
                    features.addAll(subFeats);
                }
                ii++;
            }
            return features;
        }

    }

    /**
     * Search for a feature with the given name, which has the specified aminoAcid
     * at the specified proteinPosition (1-indexed).
     *
     * @param name
     * @param proteinPosition
     * @param refAA           String symbolizing the desired amino acid
     * @param mutAA           String symbolizing the mutated amino acid
     * @return Map from genome position to features found. Feature name
     *         must be exact, but there can be multiple features with the same name
     */
    public static Map<Integer, BasicFeature> getMutation(String name, int proteinPosition, String refAA, String mutAA) {
        String nm = name.toUpperCase();
        Genome currentGenome = GENOME;
        if (!Globals.isHeadless()) {
            currentGenome = IGV.getInstance().getGenomeManager().getCurrentGenome();
        }

        Map<Integer, BasicFeature> results = new HashMap<Integer, BasicFeature>();
        List<NamedFeature> possibles = featureMap.get(nm);

        if (possibles != null) {
            synchronized (featureMap) {
                for (NamedFeature f : possibles) {
                    if (!(f instanceof BasicFeature)) {
                        continue;
                    }

                    BasicFeature bf = (BasicFeature) f;
                    Codon c = bf.getCodon(currentGenome, proteinPosition);
                    if (c == null) {
                        continue;
                    }
                    if (c.getAminoAcid().equalsByName(refAA)) {
                        Set<String> snps = AminoAcidManager.getMappingSNPs(c.getSequence(),
                                AminoAcidManager.getAminoAcidByName(mutAA));
                        if (snps.size() >= 1) {
                            results.put(c.getGenomePositions()[0], bf);
                        }
                    }
                }

            }
        }

        return results;

    }


    /**
     * Find features which have refBP as the base pair at the specified genomePosition.
     *
     * @param name
     * @param startPosition 1-based
     * @param refBP
     * @return
     */
    public static Map<Integer, BasicFeature> getMutationBP(String name, int startPosition, String refBP) {
        String nm = name.toUpperCase();
        Genome currentGenome = GENOME;
        if (!Globals.isHeadless()) {
            currentGenome = IGV.getInstance().getGenomeManager().getCurrentGenome();
        }

        Map<Integer, BasicFeature> results = new HashMap<Integer, BasicFeature>();
        List<NamedFeature> possibles = featureMap.get(nm);
        String tempBP;
        String brefBP = refBP.toUpperCase();

        if (possibles != null) {
            synchronized (featureMap) {
                for (NamedFeature f : possibles) {
                    if (!(f instanceof BasicFeature)) {
                        continue;
                    }

                    BasicFeature bf = (BasicFeature) f;
                    int genomePosition = bf.featureToGenomePosition(new int[]{startPosition - 1})[0];
                    if (genomePosition <= 0) {
                        continue;
                    }
                    final byte[] nuclSequence = currentGenome.getSequence(bf.getChr(), genomePosition, genomePosition + 1);
                    if(nuclSequence == null) {
                        continue;
                    }
                    tempBP = new String(nuclSequence);
                    if (tempBP.toUpperCase().equals(brefBP)) {
                        results.put(genomePosition, bf);
                    }
                }

            }
        }

        return results;
    }

    /**
     * Doubleton class. Can sort forward or descending, at most 2 instances.
     *
     * @param <T>
     */
    private static class FeatureComparator<T extends Feature> implements Comparator {
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

        public int compare(Object o1, Object o2) {
            T feat1 = (T) o1;
            T feat2 = (T) o2;
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
