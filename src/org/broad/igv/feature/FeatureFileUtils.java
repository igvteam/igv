/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not
 * responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL), Version 2.1 which is
 * available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */
package org.broad.igv.feature;


import org.broad.igv.Globals;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.feature.tribble.IGVBEDCodec;
import org.broad.igv.util.ParsingUtils;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;

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
     * @param iFile       - a bed file
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

                if(!chr.equals(lastChr)) {
                    lastWindowOutput = 0;
                    for(Window window : openWindows.values()) {
                        int w = window.idx;
                        int windowCenter = w * step + (windowSize / 2);
                        int windowStart = windowCenter - step/2;
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
                if(startWindow > lastWindowOutput) {
                   for(int w = lastWindowOutput; w < startWindow; w++) {
                       Window window = openWindows.get(w);
                       if (window != null) {
                           int windowCenter = w * step + (windowSize / 2);
                           int windowStart = windowCenter - step/2;
                           int windowEnd = windowStart + step;
                           pw.println(chr + "\t" + windowStart + "\t" + windowEnd + "\t" + window.count);
                       }
                   }
                    for(int w = lastWindowOutput; w < startWindow; w++) {
                        openWindows.remove(w);
                    }
                    lastWindowOutput = startWindow - 1;
                }
            }
        } finally {
            if(br != null) {
                br.close();
            }
            if(pw != null) {
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
            FeatureParser parser = AbstractFeatureParser.getInstanceFor(iFile, null);
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
            Locus locus = new Locus(tokens[1].trim());
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
