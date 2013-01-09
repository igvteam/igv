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
* Genome.java
*
* Created on November 9, 2007, 9:05 AM
*
* To change this template, choose Tools | Template Manager
* and open the template in the editor.
*/
package org.broad.igv.feature.genome;


import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.ChromosomeImpl;
import org.broad.igv.feature.Cytoband;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.ui.util.MessageUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * Simple model of a genome.  Keeps an ordered list of Chromosomes, an alias table, and genome position offsets
 * for each chromosome to support a whole-genome view.
 */
public class GenomeImpl implements Genome {

    private static Logger log = Logger.getLogger(GenomeImpl.class);
    public static final int MAX_WHOLE_GENOME = 10000;

    private String id;
    private String displayName;
    private List<String> chromosomeNames;
    private ArrayList longChromosomeNames;
    private LinkedHashMap<String, Chromosome> chromosomeMap;
    private long totalLength = -1;
    private long nominalLength = -1;
    private Map<String, Long> cumulativeOffsets = new HashMap();
    private Map<String, String> chrAliasTable;
    private Sequence sequence;
    private FeatureTrack geneTrack;

    public GenomeImpl(String id, String displayName, Sequence sequence) {
        this(id, displayName, sequence, false);
    }

    public GenomeImpl(String id, String displayName, Sequence sequence, boolean chromosOrdered) {
        this.id = id;
        this.displayName = displayName;
        this.chrAliasTable = new HashMap<String, String>();
        this.sequence = sequence;
        chromosomeNames = sequence.getChromosomeNames();

        List<Chromosome> tmpChromos = new ArrayList<Chromosome>(chromosomeNames.size());
        int maxLength = -1;
        chromosomeMap = new LinkedHashMap<String, Chromosome>(tmpChromos.size());

        for (int i = 0; i < chromosomeNames.size(); i++) {
            String chr = chromosomeNames.get(i);
            int length = sequence.getChromosomeLength(chr);
            maxLength = Math.max(maxLength, length);
            Chromosome chromo = new ChromosomeImpl(i, chr, length);
            tmpChromos.add(chromo);

            if (chromosOrdered) {
                chromosomeMap.put(chr, chromo);
            }
        }

        if (!chromosOrdered) {
            ChromosomeComparator.sortChromosomeList(tmpChromos, maxLength / 10, chromosomeMap);
            chromosomeNames = new ArrayList<String>(chromosomeMap.keySet());
        }

        initializeChromosomeAliases();
    }


    public String getChromosomeAlias(String str) {
        if (str == null) {
            return str;
        } else {
            //We intern strings used as chromosomes
            //to prevent storing multiple times
            if (!chrAliasTable.containsKey(str)) {
                chrAliasTable.put(str, str);
            }
            return chrAliasTable.get(str);
        }
    }

    public void loadUserDefinedAliases() {

        File aliasFile = new File(DirectoryManager.getGenomeCacheDirectory(), id + "_alias.tab");

        if (aliasFile.exists()) {
            if (chrAliasTable == null) chrAliasTable = new HashMap();

            BufferedReader br = null;

            try {
                br = new BufferedReader(new FileReader(aliasFile));
                addChrAliases(GenomeManager.loadChrAliases(br));
            } catch (IOException e) {
                log.error("Error loading chr alias table", e);
                if (!Globals.isHeadless())
                    MessageUtils.showMessage("<html>Error loading chromosome alias table.  Aliases will not be avaliable<br>" +
                            e.toString());
            } finally {
                if (br != null) {
                    try {
                        br.close();
                    } catch (IOException e) {
                        e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                    }
                }
            }
        }
    }

    /**
     * Pouplate the chr alias table.  The input is a collection of chromosome synonym lists.  The
     * directionality is determined by the "true" chromosome names.
     *
     * @param synonymsList
     */
    public void addChrAliases(Collection<Collection<String>> synonymsList) {

        // Convert names to a set for fast "contains" testing.
        Set<String> chrNameSet = new HashSet<String>(chromosomeNames);


        for (Collection<String> synonyms : synonymsList) {

            // Find the chromosome name as used in this genome
            String chr = null;
            for (String syn : synonyms) {
                if (chrNameSet.contains(syn)) {
                    chr = syn;
                    break;
                }
            }

            // If found register aliases
            if (chr != null) {
                for (String syn : synonyms) {
                    chrAliasTable.put(syn, chr);
                }
            } else {
                // Nothing to do.  SHould this be logged?
            }
        }
    }


    /**
     * Update the chromosome alias table with common variations
     */
    void initializeChromosomeAliases() {

        for (String name : chromosomeNames) {
            if (name.startsWith("gi|")) {
                // NCBI
                String alias = getNCBIName(name);
                chrAliasTable.put(alias, name);
            }
        }

        if (chromosomeNames.size() < 10000) {
            for (String name : chromosomeNames) {

                // UCSC Conventions
                if (name.toLowerCase().startsWith("chr")) {
                    chrAliasTable.put(name.substring(3), name);
                } else {
                    chrAliasTable.put("chr" + name, name);
                }
            }


            // These are legacy mappings,  these are now defined in the genomes alias file
            if (id.startsWith("hg") || id.equalsIgnoreCase("1kg_ref"))

            {
                chrAliasTable.put("23", "chrX");
                chrAliasTable.put("24", "chrY");
                chrAliasTable.put("MT", "chrM");
            } else if (id.startsWith("mm"))

            {
                chrAliasTable.put("21", "chrX");
                chrAliasTable.put("22", "chrY");
                chrAliasTable.put("MT", "chrM");
            } else if (id.equals("b37"))

            {
                chrAliasTable.put("chrM", "MT");
                chrAliasTable.put("chrX", "23");
                chrAliasTable.put("chrY", "24");

            }

            Collection<Map.Entry<String, String>> aliasEntries = new ArrayList(chrAliasTable.entrySet());
            for (Map.Entry<String, String> aliasEntry : aliasEntries) {
                // Illumina conventions
                String alias = aliasEntry.getKey();
                String chr = aliasEntry.getValue();
                if (!alias.endsWith(".fa")) {
                    String illuminaName = alias + ".fa";
                    chrAliasTable.put(illuminaName, chr);
                }
                if (!chr.endsWith(".fa")) {
                    String illuminaName = chr + ".fa";
                    chrAliasTable.put(illuminaName, chr);
                }
            }
        }
    }

    /**
     * Extract the user friendly name from an NCBI accession
     * example: gi|125745044|ref|NC_002229.3|  =>  NC_002229.3
     */
    public static String getNCBIName(String name) {

        String[] tokens = name.split("\\|");
        return tokens[tokens.length - 1];
    }


    public String getHomeChromosome() {
        if (getLongChromosomeNames().size() == 1 || chromosomeNames.size() > MAX_WHOLE_GENOME) {
            return getLongChromosomeNames().get(0);
        } else {
            return Globals.CHR_ALL;
        }
    }


    public Chromosome getChromosome(String chrName) {
        return chromosomeMap.get(getChromosomeAlias(chrName));
    }


    public List<String> getAllChromosomeNames() {
        return chromosomeNames;
    }


    public Collection<Chromosome> getChromosomes() {
        return chromosomeMap.values();
    }


    public long getTotalLength() {
        if (totalLength < 0) {
            totalLength = 0;
            for (Chromosome chr : chromosomeMap.values()) {
                totalLength += chr.getLength();
            }
        }
        return totalLength;
    }


    public long getCumulativeOffset(String chr) {

        Long cumOffset = cumulativeOffsets.get(chr);
        if (cumOffset == null) {
            long offset = 0;
            for (String c : getLongChromosomeNames()) {
                if (chr.equals(c)) {
                    break;
                }
                offset += getChromosome(c).getLength();
            }
            cumOffset = new Long(offset);
            cumulativeOffsets.put(chr, cumOffset);
        }
        return cumOffset.longValue();
    }

    /**
     * Covert the chromosome coordinate in BP to genome coordinates in KBP
     *
     * @param chr
     * @param locationBP
     * @return
     */
    public int getGenomeCoordinate(String chr, int locationBP) {
        return (int) ((getCumulativeOffset(chr) + locationBP) / 1000);
    }

    /**
     * Convert the genome coordinates in KBP to a chromosome coordinate
     */
    public ChromosomeCoordinate getChromosomeCoordinate(int genomeKBP) {

        long cumOffset = 0;
        List<String> wgChrNames = getLongChromosomeNames();
        for (String c : wgChrNames) {
            int chrLen = getChromosome(c).getLength();
            if ((cumOffset + chrLen) / 1000 > genomeKBP) {
                int bp = (int) (genomeKBP * 1000 - cumOffset);
                return new ChromosomeCoordinate(c, bp);
            }
            cumOffset += chrLen;
        }


        String c = wgChrNames.get(wgChrNames.size() - 1);
        int bp = (int) (genomeKBP - cumOffset) * 1000;
        return new ChromosomeCoordinate(c, bp);
    }

    /**
     * Method description
     *
     * @return
     */
    public String getId() {
        return id;
    }

    public String getNextChrName(String chr) {
        List<String> chrList = getLongChromosomeNames();
        for (int i = 0; i < chrList.size() - 1; i++) {
            if (chrList.get(i).equals(chr)) {
                return chrList.get(i + 1);
            }
        }
        return null;
    }

    public String getPrevChrName(String chr) {
        List<String> chrList = getLongChromosomeNames();
        for (int i = chrList.size() - 1; i > 0; i--) {
            if (chrList.get(i).equals(chr)) {
                return chrList.get(i - 1);
            }
        }
        return null;
    }

    public byte[] getSequence(String chr, int start, int end) {

        if (sequence == null) {
            return null;
        }

        Chromosome c = getChromosome(chr);
        if (c == null) {
            return null;
        }
        end = Math.min(end, c.getLength());
        if (end <= start) {
            return null;
        }
        return sequence.getSequence(chr, start, end);
    }

    public String getDisplayName() {
        return displayName;
    }

    public byte getReference(String chr, int pos) {
        return sequence.getBase(chr, pos);
    }


    public void setCytobands(LinkedHashMap<String, List<Cytoband>> chrCytoMap) {

        for (Map.Entry<String, List<Cytoband>> entry : chrCytoMap.entrySet()) {
            String chr = entry.getKey();
            List<Cytoband> cytobands = entry.getValue();

            Chromosome chromosome = chromosomeMap.get(chr);
            if (chromosome != null) {
                ((ChromosomeImpl) chromosome).setCytobands(cytobands);
            }
        }

    }

    public void setGeneTrack(FeatureTrack geneFeatureTrack) {
        this.geneTrack = geneFeatureTrack;
    }

    @Override
    public FeatureTrack getGeneTrack() {
        return geneTrack;
    }

    @Override
    public List<String> getLongChromosomeNames() {
        if (longChromosomeNames == null) {
            longChromosomeNames = new ArrayList(getAllChromosomeNames().size());
            long genomeLength = getTotalLength();
            for (String chrName : getAllChromosomeNames()) {
                Chromosome chr = getChromosome(chrName);
                if (chr.getLength() > (genomeLength / 3000)) {
                    longChromosomeNames.add(chrName);
                }
            }
        }
        return longChromosomeNames;
    }


    @Override
    public long getNominalLength() {
        if (nominalLength < 0) {
            nominalLength = 0;
            for (String chrName : getLongChromosomeNames()) {
                Chromosome chr = getChromosome(chrName);
                nominalLength += chr.getLength();
            }
        }
        return nominalLength;
    }

}
