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


import java.awt.Color;
import java.awt.Point;
import java.io.*;
import java.util.*;


// multiple return from loadDotBracket() and loadConnectTable()
class SeqLenAndPairs {
    public int seqLen;
    public ArrayList<Point> pairs;
    public SeqLenAndPairs(int seqLen, ArrayList<Point> pairs) {
        this.seqLen = seqLen;
        this.pairs = pairs;
    }
}

// multiple return from loadPairingProb()
class SeqLenAndBinnedPairs {
    public int seqLen;
    public  ArrayList<ArrayList<Point>> binnedPairs;
    public SeqLenAndBinnedPairs(int seqLen, ArrayList<ArrayList<Point>> binnedPairs) {
        this.seqLen = seqLen;
        this.binnedPairs = binnedPairs;
    }
}



/**
 * @author sbusan
 */
public class BasePairFileUtils {

    // TODO: support bpseq, stockholm, other formats?
    // TODO: warning dialog on file overwrite
    // FIXME: some of these might hang with empty file input

    /**
     * Convert RNA-based transcript coordinates to stranded chromosome coords for
     * base-pairing arcs.
     *
     * @param arcs
     * @param seqLen
     * @param newLeft   1-based genomic coordinate for left-most
     *                  position for input pairs sequence after transformation.
     *                  If strand is "-", this will end up being the 3-prime end
     *                  of the transcript (but the left end in genomic coords).
     * @param strand    "+" or "-"
     * @return
     */
    static LinkedList<BasePairFeature> transformArcs(LinkedList<BasePairFeature> arcs,
                                                     int seqLen,
                                                     int newLeft,
                                                     String strand) {
        LinkedList<BasePairFeature> transArcs = new LinkedList<BasePairFeature>();
        for (BasePairFeature arc : arcs) {
            String chr = arc.getChr();
            int colorIndex = arc.getColorIndex();
            int startLeft, startRight, endLeft, endRight;
            if (strand == "+") {
                startLeft = arc.getStartLeft() + newLeft - 1;
                startRight = arc.getStartRight() + newLeft - 1;
                endLeft = arc.getEndLeft() + newLeft - 1;
                endRight = arc.getEndRight() + newLeft - 1;
            } else if (strand == "-") {
                startLeft = seqLen - arc.getEndRight() + newLeft;
                startRight = seqLen - arc.getEndLeft() + newLeft;
                endLeft = seqLen - arc.getStartRight() + newLeft;
                endRight = seqLen - arc.getStartLeft() + newLeft;
            } else {
                throw new RuntimeException("Unrecognized strand (options: \"+\",\"-\")");
            }
            BasePairFeature transArc = new BasePairFeature(chr,
                    startLeft,
                    startRight,
                    endLeft,
                    endRight,
                    colorIndex);
            transArcs.add(transArc);
        }
        return transArcs;
    }

    static SeqLenAndPairs loadDotBracket(String inFile) throws
            FileNotFoundException, IOException {
        // TODO: add error messages for misformatted file
        ArrayList<Point> pairs = new ArrayList<>();
        int seqLen = 0;

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
                    int left = i+1;
                    int right = openIndices.get(k).pollLast()+1;
                    pairs.add(new Point(left, right));
                }
            }
            seqLen = struct.length();

        } finally {
            if (br != null) br.close();
        }

        return new SeqLenAndPairs(seqLen, pairs);
    }

    static SeqLenAndPairs loadConnectTable(String inFile) throws
            FileNotFoundException, IOException {
        // TODO: add error messages for misformatted file
        ArrayList<Point> pairs = new ArrayList<>();
        int seqLen = 0;

        BufferedReader br = null;

        try {
            br = new BufferedReader(new FileReader(inFile));
            seqLen = Integer.parseInt(Globals.whitespacePattern.split(br.readLine().trim())[0]);

            String nextLine;
            int n = 1;
            while ((nextLine = br.readLine()) != null && n <= seqLen) {
                String[] s = Globals.whitespacePattern.split(nextLine.trim());
                int left = Integer.parseInt(s[0]);
                int right = Integer.parseInt(s[4]);
                if (right > left) pairs.add(new Point(left, right));
                n++;
            }

        } finally {
            if (br != null) br.close();
        }

        return new SeqLenAndPairs(seqLen, pairs);
    }

    static SeqLenAndBinnedPairs loadPairingProb(String inFile) throws
            FileNotFoundException, IOException {
        // TODO: support alternate thresholds, interactive threshold update from track UI?

        ArrayList<ArrayList<Point>> binnedPairs = new ArrayList<ArrayList<Point>>();
        int seqLen = 0;

        double[] probThresh = {0.1, 0.3, 0.8};
        double[] negLogTenProbThresh = {0, 0, 0};
        for (int i = 0; i < probThresh.length; i++) {
            negLogTenProbThresh[i] = -Math.log10(probThresh[i]);
        }

        for (int i = 0; i < probThresh.length; i++) binnedPairs.add(new ArrayList<java.awt.Point>());

        BufferedReader br = null;

        try {
            br = new BufferedReader(new FileReader(inFile));
            seqLen = Integer.parseInt(Globals.whitespacePattern.split(br.readLine().trim())[0]);
            br.readLine();

            String nextLine;
            while ((nextLine = br.readLine()) != null) {
                String[] s = Globals.whitespacePattern.split(nextLine.trim());
                int left = Integer.parseInt(s[0]);
                int right = Integer.parseInt(s[1]);
                double negLogTenProb = Double.parseDouble(s[2]);
                int binIndex = -1;
                for (int i = probThresh.length - 1; i >= 0; i--) {
                    if (negLogTenProb <= negLogTenProbThresh[i]) {
                        binIndex = i;
                        break;
                    }
                }
                if (binIndex != -1) binnedPairs.get(binIndex).add(new Point(left, right));
            }

        } finally {
            if (br != null) br.close();
        }

        return new SeqLenAndBinnedPairs(seqLen, binnedPairs);
    }


    static void writeBasePairFile(String bpFile,
                                  ArrayList<Color> colors,
                                  ArrayList<String> colorLabels,
                                  ArrayList<LinkedList<BasePairFeature>> groupedArcs) throws IOException {
        PrintWriter pw = null;

        try {
            pw = new PrintWriter(new BufferedWriter(new FileWriter(bpFile)));
            // first write enumerated colors header
            for (int i=0; i<colors.size(); ++i) {
                Color color = colors.get(i);
                String label = "";
                try {
                    label = colorLabels.get(i);
                } catch (IndexOutOfBoundsException e) {
                } catch (NullPointerException e) {
                }
                pw.println("color:\t" + color.getRed() +
                        "\t" + color.getGreen() +
                        "\t" + color.getBlue() +
                        "\t" + label);
            }

            // then write arc coordinates and color index
            int colorIndex = 0;
            for (LinkedList<BasePairFeature> colorGroup : groupedArcs) {
                for (BasePairFeature arc : colorGroup) {
                    pw.println(arc.toString() + "\t" + colorIndex);
                }
                colorIndex++;
            }
        } finally {
            if (pw != null) pw.close();
        }
    }

    /**
     * Merge adjacent basepairs into helices. This makes assumptions about input pair list order.

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
                    0));
        }
        return helices;
    }

    /**
     * Convert a base pairing structure file in dot-bracket notation
     * (also known as Vienna format) to an easily parseable .bp arcs file. Does not
     * currently handle mapping coords to spliced transcripts.
     */
    public static void dotBracketToBasePairFile(String inFile,
                                                String bpFile,
                                                String chromosome,
                                                String strand,
                                                int left) throws
            FileNotFoundException, IOException {

        SeqLenAndPairs s = loadDotBracket(inFile);
        ArrayList<Point> pairs = s.pairs;
        int seqLen = s.seqLen;
        LinkedList<BasePairFeature> arcs = pairsToHelices(pairs, chromosome);
        arcs = transformArcs(arcs, seqLen, left, strand);
        ArrayList<Color> colors = new ArrayList<Color>();
        colors.add(Color.black);
        ArrayList<LinkedList<BasePairFeature>> groupedArcs = new ArrayList<LinkedList<BasePairFeature>>();
        groupedArcs.add(arcs); // list of length 1 for this case (arcs only have 1 color)
        writeBasePairFile(bpFile, colors, null, groupedArcs);
    }

    /**
     * Convert a connectivity table file as output by RNAStructure
     * to an easily parseable .bp arcs file. Does not
     * currently handle mapping coords to spliced transcripts.

     */
    public static void connectTableToBasePairFile(String inFile,
                                                  String bpFile,
                                                  String chromosome,
                                                  String strand,
                                                  int left) throws
            FileNotFoundException, IOException {

        SeqLenAndPairs s = loadConnectTable(inFile);
        ArrayList<Point> pairs = s.pairs;
        int seqLen = s.seqLen;
        LinkedList<BasePairFeature> arcs = pairsToHelices(pairs, chromosome);
        arcs = transformArcs(arcs, seqLen, left, strand);
        ArrayList<Color> colors = new ArrayList<Color>();
        colors.add(Color.black);
        ArrayList<LinkedList<BasePairFeature>> groupedArcs = new ArrayList<LinkedList<BasePairFeature>>();
        groupedArcs.add(arcs); // list of length 1 for this case (arcs only have 1 color)
        writeBasePairFile(bpFile, colors, null, groupedArcs);
    }

    /**
     * Convert a pairing probability file as output by RNAStructure
     * and/or SuperFold to an easily parseable .bp arcs file. Does not
     * currently handle mapping coords to spliced transcripts.

     */
    public static void pairingProbToBasePairFile(String inFile,
                                                 String bpFile,
                                                 String chromosome,
                                                 String strand,
                                                 int left) throws
            FileNotFoundException, IOException {

        SeqLenAndBinnedPairs s = loadPairingProb(inFile);
        ArrayList<ArrayList<Point>> binnedPairs = s.binnedPairs;
        int seqLen = s.seqLen;
        ArrayList<LinkedList<BasePairFeature>> groupedArcs = new ArrayList<LinkedList<BasePairFeature>>();
        for (ArrayList<Point> pairGroup : binnedPairs) {
            groupedArcs.add(transformArcs(pairsToHelices(pairGroup, chromosome),
                    seqLen, left, strand));
        }
        ArrayList<Color> colors = new ArrayList<Color>();
        //colors.add(new Color(255, 204, 0));
        //colors.add(new Color(72, 143, 205));
        //colors.add(new Color(81, 184, 72));
        colors.add(new Color(255, 218, 125));
        colors.add(new Color(113, 195, 209));
        colors.add(new Color(51, 114, 38));
        ArrayList<String> colorLabels = new ArrayList<String>();
        colorLabels.add("PP 10 - 30%");
        colorLabels.add("PP 30 - 80%");
        colorLabels.add("Pairing probability > 80%");
        writeBasePairFile(bpFile, colors, colorLabels, groupedArcs);
    }






}
