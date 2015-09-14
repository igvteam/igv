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


import org.broad.igv.Globals;
import org.broad.igv.feature.tribble.IGVBEDCodec;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import htsjdk.tribble.Feature;

import java.io.*;
import java.util.*;

/**
 * @author jrobinso
 */
public class FeatureFileUtils {


    /**
     * Compute feature density in units of # / window.  Note: This function could be extended to other file types by
     * using codecs.
     *
     * @param iFile      - a bed file
     * @param windowSize
     * @param step
     */
    static void computeBedDensity(String iFile, String oFile, int windowSize, int step) throws IOException {

        //FeatureCodec codec = CodecFactory.getCodec(iFile, null);

        BufferedReader br = null;
        PrintWriter pw = null;

        try {
            br = ParsingUtils.openBufferedReader(iFile);
            pw = new PrintWriter(new BufferedWriter(new FileWriter(oFile)));

            String nextLine;
            Map<Integer, Window> openWindows = new LinkedHashMap<Integer, Window>();
            int lastWindowOutput = 0;
            String lastChr = null;

            while ((nextLine = br.readLine()) != null) {
                if (nextLine.startsWith("track") || nextLine.startsWith("#")) continue;

                String[] tokens = Globals.tabPattern.split(nextLine);
                if (tokens.length < 3) continue;

                String chr = tokens[0];
                int start = Integer.parseInt(tokens[1]);
                int end = Integer.parseInt(tokens[2]);

                if (!chr.equals(lastChr)) {
                    lastWindowOutput = 0;
                    for (Window window : openWindows.values()) {
                        int w = window.idx;
                        int windowCenter = w * step + (windowSize / 2);
                        int windowStart = windowCenter - step / 2;
                        int windowEnd = windowStart + step;
                        pw.println(lastChr + "\t" + windowStart + "\t" + windowEnd + "\t" + window.count);
                    }
                    openWindows.clear();
                    lastChr = chr;
                }

                int startWindow = start / step;
                int endWindow = (end + windowSize) / step;
                for (int w = startWindow; w < endWindow; w++) {
                    Window window = openWindows.get(w);
                    if (window == null) {
                        window = new Window(w);
                        openWindows.put(w, window);
                    }
                    window.increment();
                }

                // File is sorted by start position, will never see windows < startWindow aganin
                if (startWindow > lastWindowOutput) {
                    for (int w = lastWindowOutput; w < startWindow; w++) {
                        Window window = openWindows.get(w);
                        if (window != null) {
                            int windowCenter = w * step + (windowSize / 2);
                            int windowStart = windowCenter - step / 2;
                            int windowEnd = windowStart + step;
                            pw.println(chr + "\t" + windowStart + "\t" + windowEnd + "\t" + window.count);
                        }
                    }
                    for (int w = lastWindowOutput; w < startWindow; w++) {
                        openWindows.remove(w);
                    }
                    lastWindowOutput = startWindow - 1;
                }
            }
        } finally {
            if (br != null) {
                br.close();
            }
            if (pw != null) {
                pw.close();
            }
        }

    }

    static class Window {
        int idx;
        int count;

        Window(int idx) {
            this.idx = idx;
        }

        void increment() {
            count++;
        }
    }


    public static void main(String [] args) throws IOException {
        createCanonicalGeneFile(args[0], args[1]);
    }

    /**
     * Create a "canonical" gene file from an annotation UCSC refseq type annotation file.  The "canonical" file
     * represents all the isoforms of a gene as a single feature containing the union of all exons from the
     * splice forms for the gene.  This might or might not represent an actual transcript (although it usually does).
     *
     * @param iFile
     * @param outputFile
     */
    static void createCanonicalGeneFile(String iFile, String outputFile) throws IOException {

        BufferedReader br = null;
        PrintWriter pw = null;

        try {
            br = new BufferedReader(new FileReader(iFile));
            pw = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));
            FeatureParser parser = AbstractFeatureParser.getInstanceFor(new ResourceLocator(iFile), null);
            List<Feature> features = parser.loadFeatures(br, null);
            IGVBEDCodec codec = new IGVBEDCodec();


            Map<String, List<BasicFeature>> genes = new HashMap<String, List<BasicFeature>>();
            for (Feature f : features) {

                BasicFeature transcript = (BasicFeature) f;
                String geneName = transcript.getName();

                List<BasicFeature> genelist = genes.get(geneName);
                if (genelist == null) {
                    genelist = new ArrayList<BasicFeature>();
                    genes.put(geneName, genelist);
                }

                // Loop through genes with this name to find one that overlaps
                boolean foundOverlap = false;
                for (BasicFeature gene : genelist) {
                    if (gene.overlaps(transcript)) {
                        gene.setThickStart(Math.min(gene.getThickStart(), transcript.getThickStart()));
                        gene.setThickEnd(Math.max(gene.getThickEnd(), transcript.getThickEnd()));
                        mergeExons(gene, transcript.getExons());
                        foundOverlap = true;
                        break;
                    }
                }
                if (!foundOverlap) {
                    genelist.add(transcript);

                }

            }

            for (List<BasicFeature> geneList : genes.values()) {
                for (BasicFeature gene : geneList) {
                    pw.println(codec.encode(gene));
                }
            }
        } finally {
            if (br != null) br.close();
            if (pw != null) pw.close();
        }


    }

    /**
     * Create a bed file of "TSS regions", define as the 20 bp region downstream of the start of the feature.
     * <p/>
     * It is assume that iFile is sorted by start position.
     *
     * @param iFile
     * @param outputFile
     */
    static void createTSSFile(String iFile, String outputFile) throws IOException {

        BufferedReader br = null;
        PrintWriter pw = null;

        try {
            br = new BufferedReader(new FileReader(iFile));
            pw = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));
            FeatureParser parser = AbstractFeatureParser.getInstanceFor(new ResourceLocator(iFile), null);
            List<Feature> features = parser.loadFeatures(br, null);
            IGVBEDCodec codec = new IGVBEDCodec();


            Map<String, List<BasicFeature>> genes = new HashMap<String, List<BasicFeature>>();

            int lastTSS = -1;
            for (Feature f : features) {

                BasicFeature transcript = (BasicFeature) f;

                int tss = transcript.getStrand() == Strand.POSITIVE ? f.getStart() : f.getEnd();
                if (tss != lastTSS) {
                    int tssEnd = transcript.getStrand() == Strand.POSITIVE ? tss + 20 : tss - 20;
                    pw.println(transcript.getChr() + "\t" + Math.min(tss, tssEnd) + "\t" + Math.max(tss, tssEnd));
                    lastTSS = tss;
                }


            }
        } finally {
            if (br != null) br.close();
            if (pw != null) pw.close();
        }
    }


    private static void mergeExons(BasicFeature gene, List<Exon> exons) {
        Set<IExon> exonProxies = new HashSet<IExon>(gene.getExons());
        for (Exon exon : gene.getExons()) {
            exonProxies.add(Exon.getExonProxy(exon));
        }
        for (Exon exon : exons) {
            IExon proxy = Exon.getExonProxy(exon);
            if (!exonProxies.contains(proxy)) {
                gene.addExon(exon);
            }
        }
        FeatureUtils.sortFeatureList(gene.getExons());
    }

    static void covertProbeMapToBedFile(String probeMapFile, String bedFile) throws
            FileNotFoundException, IOException {


        BufferedReader br = new BufferedReader(new FileReader(probeMapFile));
        PrintWriter pw = new PrintWriter(new FileWriter(bedFile));

        String nextLine;
        while ((nextLine = br.readLine()) != null) {
            String[] tokens = nextLine.split("\t");
            Locus locus = Locus.fromString(tokens[1].trim());
            pw.println(
                    locus.getChr() + "\t" + locus.getStart() + "\t" + locus.getEnd() + "\t" + tokens[0].trim());
        }

        br.close();
        pw.close();

    }

    static void splitEmblFileByType(String emblFile, String outputDirectory) throws IOException {

        BufferedReader br = new BufferedReader(new FileReader(emblFile));
        String nextLine;
        Set<String> codes = new HashSet();

        while ((nextLine = br.readLine()) != null) {
            if (nextLine.startsWith("FT") && (nextLine.length() > 19)) {
                String code = nextLine.substring(5, 19).trim();
                if (code.length() > 0) {
                    codes.add(code);
                }
            }
        }
        br.close();

        Map<String, PrintWriter> writers = new HashMap();
        for (String code : codes) {
            writers.put(code,
                    new PrintWriter(new FileWriter(new File(outputDirectory, code + ".embl"))));
        }

        br = new BufferedReader(new FileReader(emblFile));
        PrintWriter currentWriter = null;
        while ((nextLine = br.readLine()) != null) {
            if (nextLine.startsWith("ID")) {
                for (PrintWriter pw : writers.values()) {
                    pw.println(nextLine);
                }
            } else if (nextLine.startsWith("FT")) {
                String code = nextLine.substring(5, 19).trim();
                if (code.length() > 0) {
                    currentWriter = writers.get(code);
                }
                if (currentWriter != null) {
                    currentWriter.println(nextLine);
                }
            } else {
                currentWriter = null;
            }
        }

        br.close();
        for (PrintWriter pw : writers.values()) {
            pw.close();
        }
    }
}