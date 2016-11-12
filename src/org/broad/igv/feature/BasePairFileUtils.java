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


import org.broad.igv.feature.basepair.BasePairFeature;
import org.broad.igv.Globals;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;


import java.awt.Color;
import java.awt.Point;
import java.io.*;
import java.util.*;

/**
 * @author sbusan
 */
public class BasePairFileUtils {

    // TODO: support bpseq, stockholm, other formats?
    // TODO: warning dialog on file overwrite

    static ArrayList<Point> loadDotBracket(String inFile) throws
            FileNotFoundException, IOException {
        ArrayList<Point> pairs = new ArrayList<>();

        BufferedReader br = null;

        try {
            br = new BufferedReader(new FileReader(inFile));
            String struct = "";

            String nextLine;
            while ((nextLine = br.readLine()) != null) {

                // skip leading header lines if present
                if (nextLine.startsWith(">") || nextLine.startsWith("#")) {
                    continue;
                }

                String s = nextLine.trim();
                if (s.chars().allMatch(Character::isLetter)) {
                    // sequence line
                    continue;
                } else {
                    // assumed structure line
                    struct += s;
                }
            }

            String leftBrackets = "([{<";
            String rightBrackets = ")]}>";

            ArrayList<LinkedList<Integer>> openIndices = new ArrayList<>();
            for (int i = 0; i < leftBrackets.length(); ++i) {
                openIndices.add(new LinkedList<>());
            }

            for (int i = 0; i < struct.length(); ++i) {
                int n = leftBrackets.indexOf(struct.charAt(i));
                int k = rightBrackets.indexOf(struct.charAt(i));
                if (n >= 0) {
                    openIndices.get(n).add(i);
                } else if (k >= 0) {
                    int left = i;
                    int right = openIndices.get(k).pollLast();
                    pairs.add(new Point(left, right));
                }
            }

        } finally {
            if (br != null) br.close();
        }

        /*BufferedReader br = null;
        PrintWriter pw = null;

        try {
            br = new BufferedReader(new FileReader(iFile));
            pw = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));

            String nextLine;
            while ((nextLine = br.readLine()) != null) {

                if (nextLine.startsWith("#") || nextLine.startsWith("track") || nextLine.startsWith("browser")) {
                    pw.println(nextLine);
                    continue;
                }

                String[] tokens = Globals.whitespacePattern.split(nextLine);

                String chr = tokens[2];
                int start = Integer.parseInt(tokens[4]);
                String end = tokens[5];
                String name = tokens[12];
                String score = "1000";
                String strand = tokens[3];
                String thickStart = tokens[6];
                String thickEnd = tokens[7];
                String itemRGB = ".";
                int blockCount = Integer.parseInt(tokens[8]);

                String exonStarts = tokens[9];
                String[] stok = Globals.commaPattern.split(exonStarts);
                String[] etok = Globals.commaPattern.split(tokens[10]);

                String blockStarts = "";
                String blockSizes = "";

                for (int i = 0; i < blockCount; i++) {
                    final int bs = Integer.parseInt(stok[i]);
                    blockStarts += String.valueOf(bs - start);
                    blockSizes += String.valueOf(Integer.parseInt(etok[i]) - bs);
                    if (i != blockCount - 1) {
                        blockStarts += ",";
                        blockSizes += ",";
                    }
                }

                pw.println(chr + "\t" + start + "\t" + end + "\t" + name + "\t" + score + "\t" + strand + "\t" +
                        thickStart + "\t" + thickEnd + "\t" + itemRGB + "\t" + blockCount + "\t" + blockSizes + "\t" + blockStarts);
            }

        } finally {
            if (br != null) br.close();
            if (pw != null) pw.close();
        }*/


        return pairs;
    }

    static ArrayList<Point> loadConnectTable(String inFile) throws
            FileNotFoundException, IOException {
        ArrayList<Point> pairs = new ArrayList<>();

        return pairs;
    }

    static ArrayList<ArrayList<Point>> loadPairingProb(String inFile) throws
            FileNotFoundException, IOException {
        ArrayList<Point> pairs = new ArrayList<>();
        ArrayList<ArrayList<Point>> binnedPairs = new ArrayList<>();

        return binnedPairs;
    }


    static void writeBasePairFile(String bpFile,
                                  ArrayList<Color> colors,
                                  ArrayList<LinkedList<BasePairFeature>> groupedArcs) throws IOException {
        PrintWriter pw = null;

        try {
            pw = new PrintWriter(new BufferedWriter(new FileWriter(bpFile)));
            // first write enumerated colors header
            for (Color color : colors) {
                pw.println("color:\t" + color.getRed() + "\t" + color.getGreen() + "\t" + color.getBlue());
            }

            // then write arc coordinates
            int colorIndex = 0;
            for (LinkedList<BasePairFeature> colorGroup : groupedArcs) {
                for (BasePairFeature arc : colorGroup) {
                    pw.println(arc.toStringNoColor() + "\t" + colorIndex);
                }
                colorIndex++;
            }
        } finally {
            if (pw != null) pw.close();
        }
    }

    /**
     * Merge adjacent basepairs into helices.
     *
     * @param pairs
     * @param chromosome
     * @return
     */
    static LinkedList<BasePairFeature> pairsToHelices(ArrayList<Point> pairs,
                                                      String chromosome) {
        ArrayList<Point> bps = new ArrayList<>(pairs);
        LinkedList<LinkedList<Point>> helixPairGroups = new LinkedList<>();

        // FIXME: there should be a much faster and more elegant way to do this
        while (bps.size() > 0) {
            Point bp = bps.get(0);
            LinkedList<Point> helixPairs = new LinkedList<>();
            boolean[] removePairs = new boolean[bps.size()];
            helixPairs.add(bp);
            removePairs[0] = true;
            boolean endOfList = false;
            /* - ignore other pairs starting with this left nuc
             * - if any base pairs exist starting one nuc downstream of current left nuc
             *   and ending one nuc upstream of current right nuc, store index to remove
             *   and append to growing helix
             */
            int i = 1;
            int skippedCount = 0;
            if (i < bps.size()) {
                while (bps.get(i).x == bp.x) {
                    i++;
                    skippedCount++;
                    if (i >= bps.size()) {
                        endOfList = true;
                        break;
                    }
                }
            } else {
                endOfList = true;
            }
            while (i < bps.size()) {
                if (bps.get(i).x - bp.x > 1) {
                    // reached the end of the helix
                    helixPairGroups.add(helixPairs);
                    // remove helix pairs
                    ArrayList<Point> tmpBps = new ArrayList<>();
                    for (int k = 0; k < bps.size(); k++) {
                        if (!removePairs[k]) tmpBps.add(bps.get(k));
                    }
                    bps = tmpBps;
                    break;
                } else if (bps.get(i).y - bp.y == -1) {
                    bp = bps.get(i);
                    helixPairs.add(bp);
                    removePairs[i] = true;
                }
                i++;
                if (i >= bps.size()) {
                    endOfList = true;
                    break;
                }
            }
            if (endOfList) {
                helixPairGroups.add(helixPairs);
                // remove helix pairs
                ArrayList<Point> tmpBps = new ArrayList<>();
                for (int k = 0; k < bps.size(); k++) {
                    if (!removePairs[k]) tmpBps.add(bps.get(k));
                }
                bps = tmpBps;
            }
        }

        // convert lists of adjacent pairs to BasePairFeatures
        LinkedList<BasePairFeature> helices = new LinkedList<>();
        for (LinkedList<Point> helixPairs : helixPairGroups) {
            int startLeft = Integer.MAX_VALUE;
            int startRight = 0;
            int endLeft = Integer.MAX_VALUE;
            int endRight = 0;

            for (Point pair : helixPairs) {
                if (pair.x < startLeft) startLeft = pair.x;
                if (pair.x > startRight) startRight = pair.x;
                if (pair.y < endLeft) endLeft = pair.y;
                if (pair.y > endRight) endRight = pair.y;
            }
            helices.add(new BasePairFeature(chromosome,
                    startLeft,
                    startRight,
                    endLeft,
                    endRight,
                    null));
        }
        return helices;
    }


    public static void dotBracketToBasePairFile(String inFile,
                                                String bpFile,
                                                String chromosome,
                                                String strand,
                                                int left) throws
            FileNotFoundException, IOException {

        ArrayList<Point> pairs = loadDotBracket(inFile);
        // FIXME: transform coords here using left and strand
        LinkedList<BasePairFeature> arcs = pairsToHelices(pairs, chromosome);
        ArrayList<Color> colors = new ArrayList<>();
        colors.add(Color.black);
        ArrayList<LinkedList<BasePairFeature>> groupedArcs = new ArrayList<>();
        groupedArcs.add(arcs); // list of length 1 for this case (arcs only have 1 color)
        writeBasePairFile(bpFile, colors, groupedArcs);
    }

    public static void connectTableToBasePairFile(String inFile,
                                                  String bpFile,
                                                  String chromosome,
                                                  String strand,
                                                  int left) throws
            FileNotFoundException, IOException {

        ArrayList<Point> pairs = loadConnectTable(inFile);
    }

    public static void pairingProbToBasePairFile(String inFile,
                                                 String bpFile,
                                                 String chromosome,
                                                 String strand,
                                                 int left) throws
            FileNotFoundException, IOException {

        ArrayList<ArrayList<Point>> binnedPairs = loadPairingProb(inFile);
    }


    /**
     * Convert a base pairing structure file in dot-bracket notation
     * (also known as Vienna format) to an easily parseable .bp arcs file. Does not
     * currently handle mapping coords to spliced transcripts.
     *
     * @param dbFile        Input file
     * @param bpFile        Output file
     * @param chromosome    Associated IGV chromosome
     * @param strand        Associated strand ("+" or "-")
     * @param left          Starting left-most position (0-based)
     */


    /**
     * Convert a pairing probability file as output by RNAStructure
     * and/or SuperFold to an easily parseable .bp arcs file. Does not
     * currently handle mapping coords to spliced transcripts.
     *
     * @param dpFile        Input file
     * @param bpFile        Output file
     * @param chromosome    Associated IGV chromosome
     * @param strand        Associated strand ("+" or "-")
     * @param left          Starting left-most position (0-based)
     */

}
