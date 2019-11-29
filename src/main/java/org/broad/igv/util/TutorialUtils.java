/*
 *  The MIT License (MIT)
 *
 * Copyright (c) 2016 University of California San Diego
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 *
 */

package org.broad.igv.util;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.ChromosomeNameComparator;
import org.broad.igv.feature.genome.fasta.FastaIndexedSequence;
import org.broad.igv.sam.PicardAlignment;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;

import java.io.*;
import java.util.*;

/**
 * Utilities for creating compact files for tutorials
 */

public class TutorialUtils {


    public static void main(String[] args) throws IOException {
        sliceFasta(args[0], args[1], args[2]);
        //sampleVCF(args[0], args[1], Integer.parseInt(args[2]));
        //sliceVCF(args[0], args[1], args[2], Integer.parseInt(args[3]), Integer.parseInt(args[4]));
        //extractFasta(args[0], args[1], args[2]);
        // extractAlignments(args[0], args[1], args[2]);
        //extractFeatures(args[0], args[1], args[2]);
    }

    static void extractFasta(String inputFasta, String outputFasta, String regionsFile) throws IOException {

        FastaIndexedSequence inFasta = new FastaIndexedSequence(inputFasta);

        List<Region> regions = parseRegions(new File(regionsFile));

        PrintWriter outFasta = null;

        try {
            outFasta = new PrintWriter(new BufferedWriter(new FileWriter(outputFasta)));

            for (Region r : regions) {

                byte[] sequence = inFasta.getSequence(r.chr, r.start, r.end, true);

                outFasta.println(">" + r.name);
                outFasta.println(new String(sequence));
            }
        } finally {
            if (outFasta != null) outFasta.close();
        }
    }

    static void extractAlignments(String inputFile, String outputFile, String regionsFile)
            throws IOException {

        AlignmentReader reader = AlignmentReaderFactory.getReader(inputFile, true);
        PrintWriter out = null;

        List<Region> regions = parseRegions(new File(regionsFile));
        Map<String, List<Region>> regionMap = new HashMap<>();
        for (Region r : regions) {
            List<Region> rlist = regionMap.get(r.chr);
            if (rlist == null) {
                rlist = new ArrayList<>();
                regionMap.put(r.chr, rlist);
            }
            rlist.add(r);
        }

        try {
            out = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));
            out.println("@HD VN:1.5 SO:coordinate");
            for (Region r : regions) {
                out.println("@SQ\tSN:" + r.name + "\tLN:" + (r.end - r.start));
            }

            for (Region r : regions) {

                CloseableIterator<PicardAlignment> iter = reader.query(r.chr, r.start, r.end, false);
                while (iter.hasNext()) {
                    PicardAlignment alignment = iter.next();
                    SAMRecord record = alignment.getRecord();
                    record.setReferenceName(r.name);
                    record.setAlignmentStart(record.getAlignmentStart() - r.start);

                    if (record.getReadPairedFlag() && !record.getMateUnmappedFlag()) {
                        if (record.getMateReferenceName().equals(record.getReferenceName())) {  // todo -- handle interchr
                            record.setMateReferenceName(r.name);
                            record.setMateAlignmentStart(record.getMateAlignmentStart() - r.start);
                        } else {
                            // Try to find new mate chr
                            String newMateChr = null;
                            List<Region> rlist = regionMap.get(record.getMateReferenceName());
                            if (rlist != null) {
                                for (Region r2 : rlist) {
                                    if (r2.contains(record.getMateAlignmentStart())) {
                                        newMateChr = r.name;
                                        break;
                                    }
                                }
                            }
                            if (newMateChr != null) {
                                record.setMateReferenceName(newMateChr);
                            } else {
                                record.setMateUnmappedFlag(true);
                            }
                        }
                    }
                    out.print(record.getSAMString());
                }
                iter.close();
            }
        } finally {
            if (out != null) out.close();
        }
    }


    static void extractFeatures(String inputFile, String outputFile, String regionsFile)
            throws IOException {

        PrintWriter out = null;

        List<Region> regions = parseRegions(new File(regionsFile));

        Map<String, IntervalTree<List<Feature>>> featureMap = loadFeatures(inputFile);

        try {
            out = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));
            for (Region r : regions) {

                IntervalTree<List<Feature>> featureTree = featureMap.get(r.chr);

                if (featureTree != null) {

                    List<Interval<List<Feature>>> intervals = featureTree.findOverlapping(r.start, r.end);

                    for (Interval<List<Feature>> interval : intervals) {
                        List<Feature> features = interval.getValue();
                        for (Feature f : features) {
                            if (f.start >= r.start) {
                                String s = f.tanslate(r.name, r.start);
                                out.println(s);
                            }
                        }
                    }
                }
            }
        } finally {
            if (out != null) out.close();
        }
    }


    static List<Region> parseRegions(File file) throws IOException {

        List<Region> regions = new ArrayList<>();

        BufferedReader reader = null;
        try {
            reader = new BufferedReader(new FileReader(file));
            String nextLine;
            while ((nextLine = reader.readLine()) != null) {

                if (nextLine.startsWith("#")) continue;

                String[] tokens = Globals.whitespacePattern.split(nextLine);
                if (tokens.length > 3) {
                    regions.add(new Region(tokens[0], Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]), tokens[3]));
                }

            }
        } finally {
            if (reader != null) reader.close();
        }

        return regions;
    }

    static void sliceVCF(String inputFile, String outputFile, String chr, int start, int end) throws IOException {

        BufferedReader reader = null;
        PrintWriter out = null;


        try {
            reader =  ParsingUtils.openBufferedReader(inputFile);
            out = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));

            String nextLine;
            while ((nextLine = reader.readLine()) != null) {
                if (nextLine.startsWith("#")) {
                    out.println(nextLine);
                } else {
                    String[] tokens = Globals.tabPattern.split(nextLine);
                    String c = tokens[0];

                    if (c.equals(chr)) {
                        int pos = Integer.parseInt(tokens[1]);
                        if (pos >= start && pos <= end) {
                            out.println(nextLine);
                        }
                        if (pos > end) {
                            break;
                        }
                    }
                }
            }
        } finally {
            if (out != null) out.close();
            if (reader != null) reader.close();
        }
    }



    // Keep every nth sample (genotype).
    static void sampleVCF(String inputFile, String outputFile, int n) throws IOException {

        BufferedReader reader = null;
        PrintWriter out = null;

        try {
            reader = new BufferedReader(new FileReader(inputFile));
            out = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));

            String nextLine;
            while ((nextLine = reader.readLine()) != null) {
                if (nextLine.startsWith("##")) {
                    out.println(nextLine);
                } else {
                    String[] tokens = Globals.tabPattern.split(nextLine);
                    for(int i=0; i<9; i++) {
                        out.print(tokens[i]);
                        if(i < 8) out.print('\t');
                    }

                    for(int i=9; i<tokens.length; i+=n) {

                        if(i < tokens.length) {
                            out.print('\t');
                            out.print(tokens[i]);
                        }
                    }

                    out.println();
                }
            }
        } finally {
            if (out != null) out.close();
            if (reader != null) reader.close();
        }
    }

    static void sliceFasta(String inputFasta, String outputFasta, String chr) throws IOException {

        BufferedReader reader = null;
        PrintWriter outFasta = null;

        try {
            reader = ParsingUtils.openBufferedReader(inputFasta);
            outFasta = new PrintWriter(new BufferedWriter(new FileWriter(outputFasta)));

            String nextLine;
            boolean chrFound = false;
            while ((nextLine = reader.readLine()) != null) {

                if(chrFound) {
                    if(nextLine.startsWith(">")) {
                        break;  // Done
                    }
                    else {
                        outFasta.println(nextLine);
                    }
                }

                else {
                    if(nextLine.startsWith(">" + chr)) {
                        outFasta.println(nextLine);
                        chrFound = true;
                    }
                }
            }


        } finally {
            if (outFasta != null) outFasta.close();
        }
    }


//
//
//    private static List<Region> mergeRegions(List<Region> regions) {
//
//        List<Region> mergedRegions = new ArrayList<>(regions.size());
//
//        Region lastRegion = regions.get(0);
//        for (int i = 0; i < regions.size(); i++) {
//            Region r = regions.get(i);
//            if (r.chr.equals(lastRegion.chr)) {
//                if (r.start <= lastRegion.end) {
//                    lastRegion.end = Math.max(lastRegion.end, r.end);   // extend region
//                } else {
//                    mergedRegions.add(lastRegion);
//                }
//            } else {
//                mergedRegions.add(lastRegion);
//            }
//            lastRegion = r;
//        }
//        return mergedRegions;
//
//    }

    private static Comparator<Feature> getPositionComparator() {
        Comparator<Feature> comp = new Comparator<Feature>() {
            private Comparator<String> nameComparator = ChromosomeNameComparator.get();

            public int compare(Feature o1, Feature o2) {

                int nameComp = nameComparator.compare(o1.chr, o2.chr);
                if (nameComp == 0) {
                    return o1.start - o2.start;
                } else {
                    return nameComp;
                }
            }
        };
        return comp;
    }


    static class Region {
        String chr;
        int start;
        int end;
        String name;

        public Region(String chr, int start, int end, String name) {
            this.end = end;
            this.chr = chr;
            this.start = start;
            this.name = name;
        }

        public boolean contains(int p) {
            return p >= this.start && p <= this.end;
        }
    }

    // Hardcoded for genePred ext format (refGene.txt)
    static class Feature {

        String chr;
        int start;
        int end;
        String[] tokens;

        Feature(String[] tokens) {
            this.chr = tokens[2];
            this.start = Integer.parseInt(tokens[4]);
            this.end = Integer.parseInt(tokens[5]);
            this.tokens = tokens;
        }

        boolean overlaps(String chr, int start, int end) {
            return chr.equals(this.chr) && end > this.start && start <= this.end;
        }

        String tanslate(String newChr, int offset) {
            tokens[2] = newChr;
            tokens[4] = String.valueOf(this.start - offset);
            tokens[5] = String.valueOf(this.end - offset);
            tokens[6] = String.valueOf(Integer.parseInt(tokens[6]) - offset);
            tokens[7] = String.valueOf(Integer.parseInt(tokens[7]) - offset);

            String newExonStart = "";
            String[] exonStarts = tokens[9].split(",");
            for (String es : exonStarts) {
                newExonStart += String.valueOf(Integer.parseInt(es) - offset) + ",";
            }
            tokens[9] = newExonStart;

            String newExonEnd = "";
            String[] exonEnds = tokens[10].split(",");
            for (String es : exonEnds) {
                newExonEnd += String.valueOf(Integer.parseInt(es) - offset) + ",";
            }
            tokens[10] = newExonEnd;

            String record = tokens[0];
            for (int i = 1; i < tokens.length; i++) {
                record += "\t" + tokens[i];
            }

            return record;
        }
    }

    /*
        (
    string name;        	"Name of gene (usually transcript_id from GTF)"
    string chrom;       	"Chromosome name"
    char[1] strand;     	"+ or - for strand"
    uint txStart;       	"Transcription start position"
    uint txEnd;         	"Transcription end position"
    uint cdsStart;      	"Coding region start"
    uint cdsEnd;        	"Coding region end"
    uint exonCount;     	"Number of exons"
    uint[exonCount] exonStarts; "Exon start positions"
    uint[exonCount] exonEnds;   "Exon end positions"
    int score;            	"Score"
    string name2;       	"Alternate name (e.g. gene_id from GTF)"
    string cdsStartStat; 	"enum('none','unk','incmpl','cmpl')"
    string cdsEndStat;   	"enum('none','unk','incmpl','cmpl')"
    lstring exonFrames; 	"Exon frame offsets {0,1,2}"
    )
     */
    //585	NR_046018	chr1	+	11873	14409	14409	14409	3	11873,12612,13220,	12227,12721,14409,	0	DDX11L1	unk	unk	-1,-1,-1,
    //585	NR_024540	chr1	-	14361	29370	29370	29370	11	14361,14969,15795,16606,16857,17232,17605,17914,18267,24737,29320,	14829,15038,15947,16765,17055,17368,17742,18061,18366,24891,29370,	0	WASH7P	unk	unk	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,


    // Return a map of interval trees, keyed by chromosome name.
    static Map<String, IntervalTree<List<Feature>>> loadFeatures(String file) throws IOException {

        BufferedReader reader = null;

        reader = ParsingUtils.openBufferedReader(file);

        String nextLine;
        String lastChr = null;

        List<Feature> currentFeatureList = new ArrayList<>();
        int currentMin = Integer.MAX_VALUE;
        int currentMax = 0;

        Map<String, IntervalTree<List<Feature>>> map = new HashMap<>();

        List<Feature> features = new ArrayList<>();
        while ((nextLine = reader.readLine()) != null) {

            if (nextLine.startsWith("#") || nextLine.startsWith("track") || nextLine.startsWith("browser")) continue;

            String[] tokens = Globals.whitespacePattern.split(nextLine);
            Feature f = new Feature(tokens);
            features.add(f);
        }
        features.sort(getPositionComparator());


        for (Feature f : features) {


            if (lastChr == null) {
                currentMin = f.start;
                currentMax = f.end;
                currentFeatureList.add(f);
                IntervalTree<List<Feature>> tree = new IntervalTree<>();
                map.put(f.chr, tree);
                lastChr = f.chr;
            } else {

                if (!f.chr.equals(lastChr)) {

                    // New tree
                    IntervalTree<List<Feature>> tree = map.get(lastChr);
                    tree.insert(new Interval(currentMin, currentMax, currentFeatureList));

                    tree = new IntervalTree<>();
                    map.put(f.chr, tree);
                    lastChr = f.chr;

                    currentFeatureList = new ArrayList<>();
                    currentFeatureList.add(f);
                    currentMin = f.start;
                    currentMax = f.end;

                } else if (currentFeatureList.size() > 10) {

                    // New interval
                    IntervalTree<List<Feature>> tree = map.get(lastChr);
                    tree.insert(new Interval(currentMin, currentMax, currentFeatureList));

                    currentFeatureList = new ArrayList<>();
                    currentFeatureList.add(f);
                    currentMin = f.start;
                    currentMax = f.end;


                } else {

                    // Update interval
                    currentMin = Math.min(currentMin, f.start);
                    currentMax = Math.max(currentMax, f.end);
                    currentFeatureList.add(f);

                }
            }
        }

        return map;
    }


}
