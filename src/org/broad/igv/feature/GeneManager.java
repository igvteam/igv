/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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
 * GeneManager.java
 *
 * Created on June 12, 2007, 5:22 PM
 *
 */
package org.broad.igv.feature;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeDescriptor;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.igv.util.ResourceLocator;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

/**
 * @author jrobinso
 */
public class GeneManager {

    private static Logger log = Logger.getLogger(GeneManager.class);
    /**
     * The genome for this gene manager.  This should probably be refactored
     * so that the genome owns the gene manger, or manages the genes itself.
     */
    private Genome genome;
    Map<String, List<org.broad.tribble.Feature>> chromosomeGeneMap;
    String geneTrackName;
    Map<String, IGVFeature> geneMap;
    /**
     * Map of chromosome -> lsst gene
     */
    Map<String, Integer> maxLocationMap;
    /**
     * Map of chromosome -> longest gene
     */
    Map<String, IGVFeature> longestGeneMap;

    private TrackProperties trackProperties;

    /**
     * Creates a new instance of GeneManager
     *
     * @param genome
     */
    public GeneManager(Genome genome, String geneTrackName) {
        this.genome = genome;
        chromosomeGeneMap = new HashMap<String, List<org.broad.tribble.Feature>>();
        geneMap = new HashMap<String, IGVFeature>();
        maxLocationMap = new HashMap<String, Integer>();
        longestGeneMap = new HashMap<String, IGVFeature>();
        this.geneTrackName = geneTrackName;
    }

    public String getGeneTrackName() {
        return geneTrackName;
    }

    /**
     * Add a gene.  The gene is mapped by name and chromosome.
     *
     * @param gene
     */
    public void addGene(IGVFeature gene) {

        // If there is a genome associated with this manager only include
        // genes whose chromosome is contained in the genome.

        if (genome != null) {
            if (genome.getChromosome(gene.getChr()) == null) {
                return;
            }
        }

        // If there are multiple variant of a gene, use the longest
        IGVFeature currentGene = geneMap.get(gene.getName());
        if (gene.getIdentifier() != null) {
            geneMap.put(gene.getIdentifier().trim().toUpperCase(), gene);
        }
        if (currentGene == null) {
            if (gene.getName() != null) {
                geneMap.put(gene.getName().trim().toUpperCase(), gene);
            }
        } else {
            int w1 = currentGene.getEnd() - currentGene.getStart();
            int w2 = gene.getEnd() - gene.getStart();
            if (w2 > w1) {
                if (gene.getName() != null) {
                    geneMap.put(gene.getName().trim().toUpperCase(), gene);
                }
            }
        }


        // Update "longest gene" for this chromosome
        String chr = gene.getChr();
        List<org.broad.tribble.Feature> geneDataList = chromosomeGeneMap.get(chr);
        if (geneDataList == null) {
            geneDataList = new ArrayList<org.broad.tribble.Feature>();
            chromosomeGeneMap.put(chr, geneDataList);
            maxLocationMap.put(chr, (int) gene.getEnd());
        }
        if (gene.getEnd() > maxLocationMap.get(chr)) {
            maxLocationMap.put(chr, (int) gene.getEnd());
        }

        if (longestGeneMap.get(chr) == null) {
            longestGeneMap.put(chr, gene);
        } else {
            if (gene.getLength() > longestGeneMap.get(chr).getLength()) {
                longestGeneMap.put(chr, gene);
            }
        }

        geneDataList.add(gene);

    }

    /**
     * Method description
     *
     * @param geneName
     * @return
     */
    public IGVFeature getGene(String geneName) {
        return geneMap.get(geneName.toUpperCase());
    }

    /**
     * Return the longest gene length for the given chromosome.  Feature length is the
     * total length, including introns.  If there genes have not been loaded guess 1 MB.
     *
     * @param chr
     * @return
     */
    public int getLongestGeneLength(String chr) {
        IGVFeature longestGene = longestGeneMap.get(chr);
        return ((longestGene == null) ? 1000000 : longestGene.getLength());
    }

    /**
     * @return
     */
    public Map<String, List<org.broad.tribble.Feature>> getChromsomeGeneMap() {
        return chromosomeGeneMap;
    }

    /**
     * Return the list of genes for the specified chromosome.
     *
     * @param chromosome
     * @return
     */
    public List<org.broad.tribble.Feature> getGenesForChromosome(String chromosome) {
        List<org.broad.tribble.Feature> genes = chromosomeGeneMap.get(chromosome);
        if (genes == null) {
            log.info("No genes found for chromosome: " + chromosome);
        }
        return genes;
    }

    /**
     * Method description
     *
     * @return
     */
    public Collection<String> getChromosomes() {
        return chromosomeGeneMap.keySet();
    }

    /**
     * Return the list of genes that overlap the  specified region.
     *
     * @param chromosome
     * @param start
     * @param end
     * @return
     */
    public List<org.broad.tribble.Feature> getGenesForRegion(String chromosome, int start, int end) {
        List<org.broad.tribble.Feature> geneList = new ArrayList<org.broad.tribble.Feature>();
        for (org.broad.tribble.Feature gene : getGenesForChromosome(chromosome)) {
            if ((gene.getEnd() >= start) && (gene.getStart() <= end)) {
                geneList.add(gene);
            }
        }
        return geneList;
    }

    /**
     * Sort genes by position
     */
    public void sortGeneLists() {
        Comparator c = new Comparator() {

            public int compare(Object o1, Object o2) {
                IGVFeature g1 = (IGVFeature) o1;
                IGVFeature g2 = (IGVFeature) o2;
                return (int) (g1.getStart() - g2.getStart());
            }
        };

        for (List<org.broad.tribble.Feature> geneList : chromosomeGeneMap.values()) {
            Collections.sort(geneList, c);
        }
    }

    /**
     * Method description
     *
     * @param chr
     * @return
     */
    public int getMaximumLocation(String chr) {
        return maxLocationMap.get(chr);
    }


    /**
     * Method description
     *
     * @param genomeDescriptor
     * @return
     */
    public static AsciiLineReader getGeneReader(GenomeDescriptor genomeDescriptor) {

        InputStream is = null;
        try {

            InputStream inputStream = genomeDescriptor.getGeneStream();
            if (inputStream == null) {
                return null;
            }

            AsciiLineReader reader = null;
            if (genomeDescriptor.isGeneFileGZipFormat()) {
                is = new GZIPInputStream(inputStream);
                reader = new AsciiLineReader(is);
            } else {
                is = new BufferedInputStream(inputStream);
                reader = new AsciiLineReader(is);
            }

            return reader;
        } catch (IOException ex) {
            log.warn("Error loading the genome!", ex);
            return null;
        }

    }


    public TrackProperties getTrackProperties() {
        return trackProperties;
    }

    public void setTrackProperties(TrackProperties trackProperties) {
        this.trackProperties = trackProperties;
    }
}
