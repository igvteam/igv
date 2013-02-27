/*
 * Copyright (c) 2007-2013 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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
 *
 */
package org.broad.igv.maf;

import org.broad.igv.Globals;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.index.Interval;
import org.broad.igv.util.index.IntervalTree;
import org.broad.tribble.readers.AsciiLineReader;

import java.io.*;
import java.util.*;

/**
 * @author jrobinso
 *         Date: 2/8/13
 *         Time: 7:23 AM
 */
public class MAFIndex {

    private List<String> species;

    /**
     * Map of chromosome name -> interval tree
     */
    private Map<String, IntervalTree> intervalTrees;

    public MAFIndex() {
        intervalTrees = new HashMap<String, IntervalTree>();
    }

    /**
     * Species represented in this file.
     */
    public List<String> getSpecies() {
        return species;
    }

    public void setSpecies(List<String> species) {
        this.species = species;
    }

    /**
     * The reference species.
     */
    public String getRefId() {
        return species == null || species.isEmpty() ? null : species.get(0);
    }

    public Collection<String> getChromosomes() {
        return intervalTrees.keySet();
    }

    public IntervalTree getIntervalTree(String chr) {
        return intervalTrees.get(chr);
    }


    public void putIntervalTree(String s, IntervalTree iv) {

        intervalTrees.put(s, iv);
    }

    public void insertInterval(String lastChr, int intervalStart, int intervalEnd, long value) {
        IntervalTree iv = intervalTrees.get(lastChr);
        if (iv == null) {
            iv = new IntervalTree();
            intervalTrees.put(lastChr, iv);
        }
        iv.insert(new Interval(intervalStart, intervalEnd, value));
    }


    /**
     * @param idxFile
     * @throws java.io.IOException
     */
    public static MAFIndex loadIndex(String idxFile) throws IOException {

        MAFIndex index = new MAFIndex();
        index.species = new ArrayList<String>();

        //index = new IntervalTree<Long>();
        BufferedReader br = null;
        try {
            br = ParsingUtils.openBufferedReader(idxFile);

            // Peak at first line
            String line = br.readLine();
            if (line.startsWith("#species")) {
                // Reference species is first -- this is required.
                while ((line = br.readLine()) != null) {
                    if (line.startsWith("#end")) {
                        break;
                    }
                    index.species.add(line.trim());
                }

                IntervalTree iv = null;
                while ((line = br.readLine()) != null) {
                    if (line.trim().length() == 0) continue;
                    if (line.startsWith("#chr=")) {
                        String chr = ParsingUtils.EQ_PATTERN.split(line)[1];
                        iv = new IntervalTree();
                        index.putIntervalTree(chr, iv);
                    } else if (iv != null) {
                        String[] info = Globals.tabPattern.split(line);
                        int start = Integer.parseInt(info[0]);
                        int end = Integer.parseInt(info[1]) + start;
                        long offset = Long.parseLong(info[2]);
                        //System.out.println("Read line "+ l++);
                        iv.insert(new Interval(start, end, offset));
                    } else {
                        // log.info("Skipping line " + line);
                    }
                }

            } else {
                // A "legacy" index, load every 100 lines.
                IntervalTree iv = new IntervalTree();
                int l = 0;
                while ((line = br.readLine()) != null) {
                    if (l % 100 == 0) {
                        String[] info = Globals.tabPattern.split(line);
                        int start = Integer.parseInt(info[0]);
                        int end = Integer.parseInt(info[1]) + start;
                        long offset = Long.parseLong(info[2]);
                        //System.out.println("Read line "+ l++);
                        iv.insert(new Interval(start, end, offset));
                    }
                    l++;
                }
                index.putIntervalTree("*", iv);
            }
        } finally {
            if (br != null) br.close();
        }
        return index;
    }


    /**
     * Create an index for the MAF file.
     * <p/>
     * Example MAF lines:
     * a score=34237.000000
     * s hg19.chr1     10917 479 + 249250621 gagaggc
     * s panTro2.chr15 13606 455 - 100063
     *
     * @param alignmentFile
     * @throws IOException
     */
    public static MAFIndex createIndex(String alignmentFile) throws IOException {

        MAFIndex index = new MAFIndex();
        AsciiLineReader reader = ParsingUtils.openAsciiReader(new ResourceLocator(alignmentFile));

        String lastChr = null;

        int intervalStart = 0;
        int intervalEnd = 0;
        int blockCount = 0;
        int maxBlockCount = 50;

        String line;
        long lastOffset = 0;

        Set<String> allSpecies = new HashSet<String>();
        List<String> blockSpecies = new ArrayList<String>();
        Map<String, RunningAverage> speciesRanks = new HashMap<String, RunningAverage>();

        try {
            boolean readNext = false;
            while ((line = reader.readLine()) != null) {
                //Ignore all comment lines
                if (line.startsWith("#") || line.trim().length() == 0) {
                    continue;
                }

                if (line.startsWith("a ")) {
                    readNext = true;
                    blockCount++;
                    lastOffset = reader.getPosition();

                    // Merge species list, if any, from previous block and start new one
                    mergeSpecies(blockSpecies, allSpecies, speciesRanks);
                    blockSpecies.clear();

                } else if (line.startsWith("s ")) {

                    String[] tokens = Globals.whitespacePattern.split(line);

                    String src = tokens[1];
                    String species = src;
                    String chr = src;
                    if (src.contains(".")) {
                        String[] srcTokens = ParsingUtils.PERIOD_PATTERN.split(src);
                        species = srcTokens[0];
                        chr = srcTokens[1];
                    }
                    blockSpecies.add(species);

                    if (readNext) {
                        // This will be the reference sequence line (its always first after the "a")
                        int start = Integer.parseInt(tokens[2]);
                        int end = Integer.parseInt(tokens[3]) + start;

                        if (blockCount > maxBlockCount || (lastChr != null && !chr.equals(lastChr))) {
                            if (intervalEnd > intervalStart) {
                                index.insertInterval(lastChr, intervalStart, intervalEnd, lastOffset);
                            }
                            intervalStart = start;
                            blockCount = 1;
                        }
                        lastChr = chr;
                        intervalEnd = end;

                    }
                    readNext = false;
                } else if (line.startsWith("i ")) {
                    //We do not handle information lines yet.
                    continue;
                } else if (line.startsWith("q ")) {
                    //We do not handle quality lines yet.
                    continue;
                } else {
                    readNext = false;
                }
            }


            if (blockCount > 0) {
                index.insertInterval(lastChr, intervalStart, intervalEnd, lastOffset);
            }

            // Merge species list, if any, from previous block and start new one
            mergeSpecies(blockSpecies, allSpecies, speciesRanks);
            index.setSpecies(sortSpecies(allSpecies, speciesRanks));

            return index;


        } finally {
            reader.close();

        }

    }

    private static class RunningAverage {
        int nPts = 1;
        double average = 0;

        void addValue(double value) {
            average = (nPts * average + value) / nPts++;
        }
    }

    private static void mergeSpecies(List<String> blockSpecies, Set<String> allSpecies, Map<String, RunningAverage> speciesRank) {
        allSpecies.addAll(blockSpecies);
        for (int i = 0; i < blockSpecies.size(); i++) {
            String sp = blockSpecies.get(i);
            RunningAverage rank = speciesRank.get(sp);
            if(rank == null) {
                rank = new RunningAverage();
                speciesRank.put(sp, rank);
            }
            rank.addValue(i);
            speciesRank.put(sp, rank);
        }
        System.out.println();
    }

    private static List<String> sortSpecies(final Collection<String> allSpecies, final Map<String, RunningAverage> speciesRank) {
        List<String> speciesList = new ArrayList<String>(allSpecies);
        Collections.sort(speciesList, new Comparator<String>() {
            @Override
            public int compare(String o1, String o2) {
                double v =  speciesRank.get(o1).average - speciesRank.get(o2).average;
                return v > 0 ? 1 : (v < 0 ? -1 : 0);
            }
        });
        return speciesList;
    }


    public static void writeIndex(MAFIndex index, String indexFileName) throws IOException {

        PrintWriter pw = null;

        try {
            pw = new PrintWriter(new BufferedWriter(new FileWriter(indexFileName)));

            pw.println("#species");
            List<String> species = index.species;
            for (String sp : species) {
                pw.println(sp);
            }
            pw.println("#endSpecies");

            Collection<String> chrList = index.getChromosomes();
            for (String chr : chrList) {
                pw.println("#chr=" + chr);
                IntervalTree tree = index.getIntervalTree(chr);
                Collection<Interval> intervals = tree.getIntervals();
                for (Interval node : intervals) {
                    pw.print(String.valueOf(node.getLow()));
                    pw.print("\t");
                    pw.print(String.valueOf(node.getHigh() - node.getLow()));
                    pw.print("\t");
                    pw.println(String.valueOf(node.getValue()));
                }
            }

        } finally {
            if (pw != null) pw.close();
        }
    }

}
