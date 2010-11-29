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
import org.broad.igv.ui.util.MessageUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * Simple model of a genome
 */
public class Genome {

    private static Logger log = Logger.getLogger(Genome.class);

    private Map<String, String> chrAliasTable;

    private String id;
    private String formatVersion;
    private LinkedHashMap<String, Chromosome> chromosomeMap;
    private List<String> chromosomeNames;
    private long length = -1;
    private String annotationURL;
    private Map<String, Long> cumulativeOffsets = new HashMap();


    public Genome(String id) {
        this.id = id;
        createAliasTable();
        initAnnotationURL();
    }

    /**
     * A temporary implementation until this is supported in the genome descriptor file
     */
    static Set<String> ucscGenomes = new HashSet(Arrays.asList("hg16", "hg17", "hg18", "hg19", "mm7",
            "mm8", "mm9", "sacCer1", "sacCer2", "ce6", "canFam2", "monDom5"));

    private void initAnnotationURL() {
        if (id != null) {
            if (ucscGenomes.contains(id)) {
                annotationURL = "http://genome.ucsc.edu/cgi-bin/hgGene?hgg_gene=$$&db=" + id;
            }
        }
    }


    /**
     * Create the alias table
     * TODO -- this should be part of the genome definition file
     */
    private void createAliasTable() {

        if (chrAliasTable == null) {
            chrAliasTable = new HashMap(100);
        }

        if (id.startsWith("hg") || id.equalsIgnoreCase("1kg_ref")) {

            for (int i = 0; i < 23; i++) {
                chrAliasTable.put(String.valueOf(i), "chr" + i);
                //chrLookupTable.put("Chr" + i, "chr" + i);
            }
            chrAliasTable.put("23", "chrX");
            chrAliasTable.put("24", "chrY");
            chrAliasTable.put("chr23", "chrX");
            chrAliasTable.put("chr24", "chrY");
            chrAliasTable.put("X", "chrX");
            chrAliasTable.put("Y", "chrY");
            chrAliasTable.put("M", "chrM");
            chrAliasTable.put("MT", "chrM");
        } else if (id.startsWith("mm")) {

            for (int i = 0; i < 21; i++) {
                chrAliasTable.put(String.valueOf(i), "chr" + i);
                //chrLookupTable.put("Chr" + i, "chr" + i);
            }
            chrAliasTable.put("21", "chrX");
            chrAliasTable.put("22", "chrY");
            chrAliasTable.put("chr21", "chrX");
            chrAliasTable.put("chr22", "chrY");
            chrAliasTable.put("X", "chrX");
            chrAliasTable.put("Y", "chrY");
            chrAliasTable.put("M", "chrM");
            chrAliasTable.put("MT", "chrM");
        } else if (id.equals("b37")) {
            for (int i = 0; i < 23; i++) {
                chrAliasTable.put("chr" + i, String.valueOf(i));
            }
            chrAliasTable.put("chrX", "X");
            chrAliasTable.put("chrY", "Y");
            chrAliasTable.put("chrM", "MT");

        }

    }


    public String getChromosomeAlias(String str) {
        if (chrAliasTable == null) {
            return str;
        }

        String chr = chrAliasTable.get(str);
        if (chr == null) {
            if (ucscGenomes.contains(id) && !str.startsWith("chr") || str.equals(Globals.CHR_ALL)) {
                chr = "chr" + str;
            } else {
                chr = str;
            }

            chrAliasTable.put(str, chr);
        }
        return chr;
    }

    public void loadUserDefinedAliases() {

        File aliasFile = new File(Globals.getGenomeCacheDirectory(), id + "_alias.tab");

        if (aliasFile.exists()) {
            if(chrAliasTable == null) chrAliasTable = new HashMap();

            BufferedReader br = null;

            try {
                br = new BufferedReader(new FileReader(aliasFile));
                loadChrAliases(br);
            } catch (Exception e) {
                log.error("Error loading chr alias table", e);
                MessageUtils.showMessage("<html>Error loading chromosome alias table.  Aliases will not be avaliable<br>" +
                        e.toString());
            }

            finally {
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

    public void loadChrAliases(BufferedReader br) throws IOException {
        String nextLine = "";
        while ((nextLine = br.readLine()) != null) {
            String[] kv = nextLine.split("\t");
            if (kv.length > 1) {
                chrAliasTable.put(kv[0], kv[1]);
            }
        }
    }


    public void setChromosomeMap(LinkedHashMap<String, Chromosome> chromosomeMap, boolean chromosomesAreOrdered) {
        this.chromosomeMap = chromosomeMap;
        this.chromosomeNames = new LinkedList<String>(chromosomeMap.keySet());
        if (!chromosomesAreOrdered) {
            Collections.sort(chromosomeNames, new ChromosomeComparator());
        }
    }


    public String getHomeChromosome() {
        if (getChromosomeNames().size() == 1 || chromosomeNames.size() > 5000) {
            return getChromosomeNames().get(0);
        } else {
            return Globals.CHR_ALL;
        }
    }


    public Chromosome getChromosome(String chrName) {
        return chromosomeMap.get(chrName);
    }


    public List<String> getChromosomeNames() {
        return chromosomeNames;
    }


    public Collection<Chromosome> getChromosomes() {
        return chromosomeMap.values();
    }


    public long getLength() {
        if (length < 0) {
            length = 0;
            for (Chromosome chr : chromosomeMap.values()) {
                length += chr.getLength();
            }
        }
        return length;
    }


    public long getCumulativeOffset(String chr) {

        Long cumOffset = cumulativeOffsets.get(chr);
        if (cumOffset == null) {
            long offset = 0;
            for (String c : getChromosomeNames()) {
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
        for (String c : chromosomeNames) {
            int chrLen = getChromosome(c).getLength();
            if ((cumOffset + chrLen) / 1000 > genomeKBP) {
                int bp = (int) (genomeKBP * 1000 - cumOffset);
                return new ChromosomeCoordinate(c, bp);
            }
            cumOffset += chrLen;
        }

        String c = chromosomeNames.get(chromosomeNames.size() - 1);
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


    public String getAnnotationURL() {
        return annotationURL;
    }

    public String getNextChrName(String chr) {
        List<String> chrList = getChromosomeNames();
        for (int i = 0; i < chrList.size() - 1; i++) {
            if (chrList.get(i).equals(chr)) {
                return chrList.get(i + 1);
            }
        }
        return null;
    }

    public String getPrevChrName(String chr) {
        List<String> chrList = getChromosomeNames();
        for (int i = chrList.size() - 1; i > 0; i--) {
            if (chrList.get(i).equals(chr)) {
                return chrList.get(i - 1);
            }
        }
        return null;
    }


    public static class ChromosomeCoordinate {
        private String chr;
        private int coordinate;

        public ChromosomeCoordinate(String chr, int coordinate) {
            this.chr = chr;
            this.coordinate = coordinate;
        }

        public String getChr() {
            return chr;
        }

        public int getCoordinate() {
            return coordinate;
        }
    }

    /**
     * Comparator for chromosome names.
     */
    public static class ChromosomeComparator implements Comparator<String> {

        /**
         * @param chr1
         * @param chr2
         * @return
         */
        public int compare(String chr1, String chr2) {

            try {

                // Special rule -- put the mitochondria at the end
                if (chr1.equals("chrM") || chr1.equals("MT")) {
                    return 1;
                } else if (chr2.equals("chrM") || chr2.equals("MT")) {
                    return -1;
                }

                // Find the first digit
                int idx1 = findDigitIndex(chr1);
                int idx2 = findDigitIndex(chr2);
                if (idx1 == idx2) {
                    String alpha1 = idx1 == -1 ? chr1 : chr1.substring(0, idx1);
                    String alpha2 = idx2 == -1 ? chr2 : chr2.substring(0, idx2);
                    int alphaCmp = alpha1.compareTo(alpha2);
                    if (alphaCmp != 0) {
                        return alphaCmp;
                    } else {
                        int dig1 = Integer.parseInt(chr1.substring(idx1));
                        int dig2 = Integer.parseInt(chr2.substring(idx2));
                        return dig1 - dig2;
                    }
                } else if (idx1 == -1) {
                    return 1;

                } else if (idx2 == -1) {
                    return -1;
                }
                return idx1 - idx2;
            }
            catch (Exception numberFormatException) {
                return 0;
            }

        }

        int findDigitIndex(String chr) {

            int n = chr.length() - 1;
            if (!Character.isDigit(chr.charAt(n))) {
                return -1;
            }

            for (int i = n - 1; i > 0; i--) {
                if (!Character.isDigit(chr.charAt(i))) {
                    return i + 1;
                }
            }
            return 0;
        }


    }
}
