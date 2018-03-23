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

package org.broad.igv.maf;

import org.broad.igv.Globals;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.index.Interval;
import org.broad.igv.util.index.IntervalTree;
import htsjdk.tribble.readers.AsciiLineReader;

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

    /**
     * The # of alignments represented by an interval in the tree.
     * Note: This is not private so that it can be manipulated by unit tests.
     */
    public static int blockSize = 50;

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

        IntervalTree iv = intervalTrees.get(chr);
        if(iv == null) {
            iv = intervalTrees.get("*"); // To support legacy MAF indeces, files are split by chromosome
        }
        return iv;
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
                        iv.insert(new Interval(start, end, offset));
                    } else {
                        // log.info("Skipping line " + line);
                    }
                }

            } else {
                // A "legacy" index, created for Broad hosted files that are separated by chromosome.
                // Every alignment is indexed, which is overkill.  Below we lump them into blocks of 50.
                IntervalTree iv = new IntervalTree();
                int l = 0;
                int intervalStart = 0;
                int intervalEnd = 0;
                long lastOffset = 0;
                while ((line = br.readLine()) != null) {
                    String[] info = Globals.tabPattern.split(line);
                    int start = Integer.parseInt(info[0]);
                    intervalEnd = Integer.parseInt(info[1]) + start;
                    if (l % 50 == 0) {
                        iv.insert(new Interval(intervalStart, intervalEnd, lastOffset));
                        intervalStart = intervalEnd;
                        lastOffset = Long.parseLong(info[2]);
                    }
                    l++;
                }

                if(intervalEnd > intervalStart) {
                    iv.insert(new Interval(intervalStart, intervalEnd, lastOffset));
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


        try {

            String lastChr = null;
            int intervalStart = 0;
            int intervalEnd = 0;
            int blockCount = 0;
            String line;
            long lastOffset = 0;
            boolean newBlock = false;

            Set<String> allSpecies = new HashSet<String>();
            List<String> blockSpecies = new ArrayList<String>();
            Map<String, RunningAverage> speciesRanks = new HashMap<String, RunningAverage>();

            while ((line = reader.readLine()) != null) {
                //Ignore all comment lines
                if (line.startsWith("#") || line.trim().length() == 0) {
                    continue;
                }

                if (line.startsWith("a ")) {
                    newBlock = true;
                    blockCount++;

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
                    if (lastChr == null) lastChr = chr;

                    blockSpecies.add(species);

                    if (newBlock) {
                        // This will be the reference sequence line (its always first after the "a")
                        int start = Integer.parseInt(tokens[2]);
                        int end = Integer.parseInt(tokens[3]) + start;

                        if ((!chr.equals(lastChr) && blockCount > 0) ||
                                blockCount > blockSize) {

                            // Record previous interval and start a new one.
                            index.insertInterval(lastChr, intervalStart, intervalEnd, lastOffset);

                            blockCount = 1;
                            lastOffset = reader.getPosition();
                            intervalStart = start;
                        }

                        lastChr = chr;
                        intervalEnd = end;
                        newBlock = false;
                    }

                } else if (line.startsWith("i ")) {
                    //We do not handle information lines yet.
                    continue;
                } else if (line.startsWith("q ")) {
                    //We do not handle quality lines yet.
                    continue;
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
            if (rank == null) {
                rank = new RunningAverage();
                speciesRank.put(sp, rank);
            }
            rank.addValue(i);
            speciesRank.put(sp, rank);
        }
    }

    private static List<String> sortSpecies(final Collection<String> allSpecies, final Map<String, RunningAverage> speciesRank) {
        List<String> speciesList = new ArrayList<String>(allSpecies);
        Collections.sort(speciesList, new Comparator<String>() {
            @Override
            public int compare(String o1, String o2) {
                double v = speciesRank.get(o1).average - speciesRank.get(o2).average;
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
