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

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ui.IGV;

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
    private static Map<String, NamedFeature> featureMap = Collections.synchronizedSortedMap(new TreeMap<String, NamedFeature>());
    private static Genome GENOME = null;

    public static void addFeature(NamedFeature feature) {

        final String name = feature.getName();
        if (name != null && name.length() > 0) {
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

    public static void put(String name, NamedFeature feature) {

        String key = name.toUpperCase();

        Genome currentGenome = GENOME;
        if (!Globals.isHeadless()) {
            currentGenome = IGV.getInstance().getGenomeManager().getCurrentGenome();
        }

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
    }

    static void setGenome(Genome genome) {
        GENOME = genome;
    }


    public static void addFeature(String name, NamedFeature feature) {
        if (Globals.isHeadless()) {
            return;
        }
        featureMap.put(name.toUpperCase(), feature);

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

    /**
     * Return the feature, if any, with the given name.  Genes are given
     * precedence.
     */
    public static NamedFeature getFeature(String nm) {
        String name = nm.trim().toUpperCase();
        return featureMap.get(name);
    }

    /**
     * Get all features which match nm. Not necessarily
     * an exact match. Currenty implementation will match anything
     * for which nm is at the beginning, including but not limited to
     * exact matches.
     * <p/>
     * NOTE: "It is imperative that the user manually synchronize
     * on [this sorted map] when iterating over any of its
     * collection views, or the collections views of any of its
     * subMap, headMap or tailMap views". See
     * <a href="http://docs.oracle.com/javase/6/docs/api/java/util/Collections.html#synchronizedSortedMap%28java.util.SortedMap%29> here</a>
     *
     * @param nm : Search string. Features which begin with this
     *           string will be found.
     * @return
     */
    public static Map<String, NamedFeature> getFeatures(String nm) {
        String name = nm.trim().toUpperCase();
        SortedMap<String, NamedFeature> treeMap = (SortedMap) featureMap;
        //Search is inclusive to first argument, exclusive to second
        return treeMap.subMap(name, name + (char) ((int) 'Z' + 1));
    }

}
