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
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.util.collections.MultiMap;
import htsjdk.tribble.Feature;

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
    private static final int MAX_DUPLICATE_COUNT = 20;

    public static void addFeature(NamedFeature feature, Genome genome) {

        final String name = feature.getName();
        if (name != null && name.length() > 0 && !name.equals(".")) {
            put(name, feature, genome);
        }
        if (feature instanceof IGVFeature) {
            final IGVFeature igvFeature = (IGVFeature) feature;
            final String id = igvFeature.getIdentifier();
            if (id != null && id.length() > 0) {
                put(id, feature, genome);
            }

            addByAttributes(igvFeature, genome);

            List<Exon> exons = igvFeature.getExons();
            if (exons != null) {
                for (Exon exon : exons) {
                    addByAttributes(exon, genome);
                }
            }
        }
    }

    private static void addByAttributes(IGVFeature igvFeature, Genome genome) {
        MultiMap<String, String> attributes = igvFeature.getAttributes();
        if (attributes != null) {
            for (String value : attributes.values()) {
                if (value.length() < 20) {
                    put(value, igvFeature, genome);
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
     * @param genome The genome which these features belong to. Used for checking chromosomes
     * @return true if successfully added, false if not
     */
    static boolean put(String name, NamedFeature feature, Genome genome) {
        String key = name.toUpperCase();
        if (!Globals.isHeadless()) {
            Genome currentGenome = genome != null ? genome : GenomeManager.getInstance().getCurrentGenome();
            if (currentGenome != null && currentGenome.getChromosome(feature.getChr()) == null) {
                return false;
            }
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

                return currentList.add(feature);

            }
        }
    }

    /*
        String key = name.toUpperCase();

        Genome currentGenome = GenomeManager.getInstance().getCurrentGenome();
        if (currentGenome == null || currentGenome.getChromosome(feature.getChr()) != null) {
            NamedFeature currentFeature = featureMap.get(key);
            if (currentFeature == null) {
                featureMap.put(key, feature);
            } else {
                // If there are multiple features, prefer the one that is NOT on a "random" chromosome.
                // This is a hack, but an important one for the human assemblies
                String featureChr = feature.getChr().toLowerCase();
                String currentFeatureChr = currentFeature.getChr();
                if (featureChr.contains("random") || featureChr.contains("chrun") || featureChr.contains("hap")) {
                    return;
                } else if (currentFeatureChr.contains("random") || currentFeatureChr.contains("chrun") ||
                        currentFeatureChr.contains("hap")) {
                    featureMap.put(key, feature);
                    return;
                }

                // If there are multiple features, use or keep the longest one
                int w1 = currentFeature.getEnd() - currentFeature.getStart();
                int w2 = feature.getEnd() - feature.getStart();
                if (w2 > w1) {
                    featureMap.put(key, feature);
                }

            }

        }

     */


    public static void addFeature(String name, NamedFeature feature, Genome genome) {
        put(name.toUpperCase(), feature, genome);
    }


    private FeatureDB() {
    }


    public static void addFeatures(List<htsjdk.tribble.Feature> features, Genome genome) {
        for (htsjdk.tribble.Feature feature : features) {
            if (feature instanceof IGVFeature)
                addFeature((IGVFeature) feature, genome);
        }
    }


    public static void clearFeatures() {
        featureMap.clear();
    }

    static int size() {
        return featureMap.size();
    }

    /**
     * Return a feature with the given name.
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

    /**
     * Shortcut to getFeaturesList(name, limit, true)
     *
     * @param name
     * @param limit
     * @return
     * @see #getFeaturesList(String, int, boolean)
     */
    public static List<NamedFeature> getFeaturesList(String name, int limit) {
        return getFeaturesList(name, limit, true);
    }

    /**
     * Get a list of features which start with the provided name.
     * Note that matches can be inexact
     *
     * @param name
     * @param limit
     * @param longestOnly Whether to take only the longest feature for each name
     * @return
     */
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
     * at the specified (1-indexed) proteinPosition .
     *
     * @param name
     * @param proteinPosition 1-indexed position along the feature in which the amino acid is found, in protein coordinates
     * @param refAA           String symbolizing the desired amino acid
     * @param mutAA           String symbolizing the mutated amino acid
     * @param currentGenome
     * @return Map from genome position to features found. Feature name
     *         must be exact, but there can be multiple features with the same name
     */
    public static Map<Integer, BasicFeature> getMutationAA(String name, int proteinPosition, String refAA,
                                                           String mutAA, Genome currentGenome) {
        String nm = name.toUpperCase();

        if (!Globals.isHeadless() && currentGenome == null) {
            currentGenome = GenomeManager.getInstance().getCurrentGenome();
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
                        Set<String> snps = AminoAcidManager.getInstance().getMappingSNPs(c.getSequence(),
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
     * Find features with a given name, which have refNT as the base pair at the specified position within the feature.
     * refNT considered based on the read strand, so a negative strand feature with A at position 1 on the positive strand
     * would be found only if refNT = T.
     *
     * @param name          Feature name
     * @param startPosition 1-based location within the feature
     * @param refNT         Nucleotide (A, G, C, T) of feature.
     * @param currentGenome The genome in which to search
     * @return
     */
    public static Map<Integer, BasicFeature> getMutationNT(String name, int startPosition, String refNT, Genome currentGenome) {
        String nm = name.toUpperCase();
        if (!Globals.isHeadless() && currentGenome == null) {
            currentGenome = GenomeManager.getInstance().getCurrentGenome();
        }

        Map<Integer, BasicFeature> results = new HashMap<Integer, BasicFeature>();
        List<NamedFeature> possibles = featureMap.get(nm);
        String tempNT;
        String brefNT = refNT.toUpperCase();

        if (possibles != null) {
            synchronized (featureMap) {
                for (NamedFeature f : possibles) {
                    if (!(f instanceof BasicFeature)) {
                        continue;
                    }

                    BasicFeature bf = (BasicFeature) f;

                    int genomePosition = bf.featureToGenomePosition(new int[]{startPosition - 1})[0];
                    if (genomePosition < 0) {
                        continue;
                    }
                    final byte[] nuclSequence = currentGenome.getSequence(bf.getChr(), genomePosition, genomePosition + 1);
                    if (nuclSequence == null) {
                        continue;
                    }
                    tempNT = new String(nuclSequence);
                    if (bf.getStrand() == Strand.NEGATIVE) {
                        tempNT = AminoAcidManager.getNucleotideComplement(tempNT);
                    }

                    if (tempNT.toUpperCase().equals(brefNT)) {
                        results.put(genomePosition, bf);
                    }
                }

            }
        }

        return results;
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
