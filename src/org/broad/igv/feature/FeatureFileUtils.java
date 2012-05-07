/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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
package org.broad.igv.feature;

//~--- JDK imports ------------------------------------------------------------

import org.broad.igv.feature.tribble.IGVBEDCodec;
import org.broad.tribble.Feature;

import java.io.*;
import java.util.*;

/**
 * chromosome
 * coding_exon
 * contig
 * EST
 * EST_match
 * exon
 * expressed_sequence_match
 * gap
 * gene
 * ##gff-version 3
 * match
 * microarray_oligo
 * mRNA
 * processed_transcript
 * region
 *
 * @author jrobinso
 */
public class FeatureFileUtils {

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
