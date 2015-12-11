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
* Genome.java
*
* Created on November 9, 2007, 9:05 AM
*
* To change this template, choose Tools | Template Manager
* and open the template in the editor.
*/
package org.broad.igv.feature.genome;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.Cytoband;
import org.broad.igv.track.FeatureTrack;

import java.io.*;
import java.util.*;

/**
 * Simple model of a genome.  Keeps an ordered list of Chromosomes, an alias table, and genome position offsets
 * for each chromosome to support a whole-genome view.
 */
public class Genome {

    private static Logger log = Logger.getLogger(Genome.class);
    public static final int MAX_WHOLE_GENOME = 10000;

    private String id;
    private String displayName;
    private List<String> chromosomeNames;
    private ArrayList<String> longChromosomeNames;
    private LinkedHashMap<String, Chromosome> chromosomeMap;
    private long totalLength = -1;
    private long nominalLength = -1;
    private Map<String, Long> cumulativeOffsets = new HashMap();
    private Map<String, String> chrAliasTable;
    private Sequence sequence;
    private FeatureTrack geneTrack;
    private String species;

    /**
     * @param id
     * @param displayName
     * @param sequence       the reference Sequence object.  Can be null.
     * @param chromosOrdered Whether the chromosomes are already ordered. If false, they will be sorted.
     */
    public Genome(String id, String displayName, Sequence sequence, boolean chromosOrdered) {
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
            Chromosome chromo = new Chromosome(i, chr, length);
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


    /**
     * Alternate constructor for defining a minimal genome, usually from parsing a chrom.sizes file.
     *
     * @param id
     * @param chromosomes
     */
    public Genome(String id, List<Chromosome> chromosomes) {
        this.id = id;
        this.displayName = id;
        this.chrAliasTable = new HashMap<String, String>();
        this.sequence = null;

        chromosomeNames = new ArrayList<String>(chromosomes.size());
        chromosomeMap = new LinkedHashMap<String, Chromosome>(chromosomes.size());
        for (Chromosome chromosome : chromosomes) {
            chromosomeNames.add(chromosome.getName());
            chromosomeMap.put(chromosome.getName(), chromosome);
        }
        initializeChromosomeAliases();

    }


    public String getCanonicalChrName(String str) {
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

    public Map<String, String> getChrAliasTable() {
        return chrAliasTable;
    }


    /**
     * Populate the chr alias table.  The input is a collection of chromosome synonym lists.  The
     * directionality is determined by the "true" chromosome names.
     *
     * @param synonymsList
     */
    public void addChrAliases(Collection<Collection<String>> synonymsList) {

        if(chrAliasTable == null) chrAliasTable = new HashMap<String, String>();

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
        chrAliasTable.putAll(getAutoAliases());

    }


    Map<String, String> getAutoAliases() {

        Map<String, String> autoAliases = new HashMap<String, String>();

        for (String name : chromosomeNames) {
            if (name.startsWith("gi|")) {
                // NCBI
                String alias = getNCBIName(name);
                autoAliases.put(alias, name);

                // Also strip version number out, if present
                int dotIndex = alias.lastIndexOf('.');
                if(dotIndex > 0) {
                    alias = alias.substring(0, dotIndex);
                    autoAliases.put(alias, name);
                }
            }
        }

        if (chromosomeNames.size() < 10000) {
            for (String name : chromosomeNames) {

                // UCSC Conventions
                if (name.toLowerCase().startsWith("chr")) {
                    autoAliases.put(name.substring(3), name);
                } else {
                    autoAliases.put("chr" + name, name);
                }
            }


            // These are legacy mappings,  these are now defined in the genomes alias file
            if (id.startsWith("hg") || id.equalsIgnoreCase("1kg_ref"))

            {
                autoAliases.put("23", "chrX");
                autoAliases.put("24", "chrY");
                autoAliases.put("MT", "chrM");
            } else if (id.startsWith("mm"))

            {
                autoAliases.put("21", "chrX");
                autoAliases.put("22", "chrY");
                autoAliases.put("MT", "chrM");
            } else if (id.equals("b37"))

            {
                autoAliases.put("chrM", "MT");
                autoAliases.put("chrX", "23");
                autoAliases.put("chrY", "24");

            }

            Collection<Map.Entry<String, String>> aliasEntries = new ArrayList(autoAliases.entrySet());
            for (Map.Entry<String, String> aliasEntry : aliasEntries) {
                // Illumina conventions
                String alias = aliasEntry.getKey();
                String chr = aliasEntry.getValue();
                if (!alias.endsWith(".fa")) {
                    String illuminaName = alias + ".fa";
                    autoAliases.put(illuminaName, chr);
                }
                if (!chr.endsWith(".fa")) {
                    String illuminaName = chr + ".fa";
                    autoAliases.put(illuminaName, chr);
                }
            }
        }
        return autoAliases;
    }
    /**
     * Extract the user friendly name from an NCBI accession
     * example: gi|125745044|ref|NC_002229.3|  =>  NC_002229.3
     */
    public static String getNCBIName(String name) {

        String[] tokens = name.split("\\|");
        return tokens[tokens.length - 1];
    }


    /**
     * Return the chromosome name associated with the "home" button,  usually the whole genome chromosome.
     *
     * @return
     */
    public String getHomeChromosome() {
        if (chromosomeNames.size() == 1 || chromosomeNames.size() > MAX_WHOLE_GENOME) {
            return chromosomeNames.get(0);
        } else {
            return Globals.CHR_ALL;
        }
    }


    public Chromosome getChromosome(String chrName) {
        return chromosomeMap.get(getCanonicalChrName(chrName));
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
     * Covert the chromosome coordinate in basepairs to genome coordinates in kilo-basepairs
     *
     * @param chr
     * @param locationBP
     * @return The overall genome coordinate, in kilo-bp
     */
    public int getGenomeCoordinate(String chr, int locationBP) {
        return (int) ((getCumulativeOffset(chr) + locationBP) / 1000);
    }

    /**
     * Translate a genome coordinate, in kilo-basepairs, to a chromosome & position in basepairs.
     *
     * @param genomeKBP The "genome coordinate" in kilo-basepairs.  This is the distance in kbp from the start of the
     *                  first chromosome.
     * @return the position on the corresponding chromosome
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

    public String getSpecies() {
        if (species == null) {
            species = Genome.getSpeciesForID(id);
        }
        return species;
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

    /**
     * Return the nucleotide sequence on the + strand for the genomic interval.  This method can return null
     * if sequence is not available.
     *
     * @param chr
     * @param start start position in "zero-based" coordinates
     * @param end   end position
     * @return sequence, or null if not available
     * @api
     */
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

    /**
     * Return the reference base at the given position.  Can return null if reference sequence is unknown
     *
     * @param chr
     * @param pos
     * @return the reference base, or null if unknown
     */
    public byte getReference(String chr, int pos) {
        return sequence == null ? null : sequence.getBase(chr, pos);
    }


    public void setCytobands(LinkedHashMap<String, List<Cytoband>> chrCytoMap) {

        for (Map.Entry<String, List<Cytoband>> entry : chrCytoMap.entrySet()) {
            String chr = entry.getKey();
            List<Cytoband> cytobands = entry.getValue();

            Chromosome chromosome = chromosomeMap.get(chr);
            if (chromosome != null) {
                chromosome.setCytobands(cytobands);
            }
        }

    }

    public void setGeneTrack(FeatureTrack geneFeatureTrack) {
        this.geneTrack = geneFeatureTrack;
    }

    /**
     * Return the annotation track associated with this genome.   Can return null
     *
     * @return a FeatureTrack, or null
     */
    public FeatureTrack getGeneTrack() {
        return geneTrack;
    }

    /**
     * Return "getChromosomeNames()" with small chromosomes removed.
     *
     * @return
     */
    public List<String> getLongChromosomeNames() {
        if (longChromosomeNames == null) {
            longChromosomeNames = new ArrayList<String>(getAllChromosomeNames().size());
            long genomeLength = getTotalLength();
            int maxChromoLength = -1;
            for (String chrName : getAllChromosomeNames()) {
                Chromosome chr = getChromosome(chrName);
                int length = chr.getLength();
                maxChromoLength = Math.max(maxChromoLength, length);
                if (length > (genomeLength / 3000)) {
                    longChromosomeNames.add(chrName);
                }
            }

            /**
             * At this point, we should have some long chromosome names.
             * However, some genomes (draft ones perhaps) maybe have many small ones
             * which aren't big enough. We arbitrarily take those which are above
             * half the size of the max, only if the first method didn't work.
             */
            if (longChromosomeNames.size() == 0) {
                for (String chrName : getAllChromosomeNames()) {
                    Chromosome chr = getChromosome(chrName);
                    int length = chr.getLength();
                    if (length > maxChromoLength / 2) {
                        longChromosomeNames.add(chrName);
                    }
                }
            }
        }
        return longChromosomeNames;
    }

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


    // TODO A hack (
    // obviously),  we need to record a species in the genome definitions
    private static Map<String, String> ucscSpeciesMap;

    private static synchronized String getSpeciesForID(String id) {
        if (ucscSpeciesMap == null) {
            ucscSpeciesMap = new HashMap<String, String>();

            InputStream is = null;

            try {
                is = Genome.class.getResourceAsStream("speciesMapping.txt");
                BufferedReader br = new BufferedReader(new InputStreamReader(is));

                String nextLine;
                while ((nextLine = br.readLine()) != null) {
                    if (nextLine.startsWith("#")) continue;
                    String[] tokens = Globals.whitespacePattern.split(nextLine);
                    if (tokens.length == 2) {
                        ucscSpeciesMap.put(tokens[0], tokens[1]);
                    } else {
                        log.error("Unexpected number of tokens in species mapping file for line: " + nextLine);
                    }
                }
            } catch (IOException e) {
                log.error("Error reading species mapping table", e);
            } finally {
                if (is != null) try {
                    is.close();
                } catch (IOException e) {
                    log.error("", e);
                }
            }

        }

        for (Map.Entry<String, String> entry : ucscSpeciesMap.entrySet()) {
            if (id.startsWith(entry.getKey())) {
                return entry.getValue();
            }
        }
        return null;
    }


}
