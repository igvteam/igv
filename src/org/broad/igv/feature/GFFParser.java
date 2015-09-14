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

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.exceptions.ParserException;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.GFFCodec;
import org.broad.igv.renderer.GeneTrackRenderer;
import org.broad.igv.renderer.IGVFeatureRenderer;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;

import java.io.*;
import java.util.*;

/**
 * @deprecated use org.broad.igv.track.GFFFeatureSource
 * User: jrobinso
 */


public class GFFParser implements FeatureParser {

    static Logger log = Logger.getLogger(GFFParser.class);

    private TrackProperties trackProperties = null;

    public List<htsjdk.tribble.Feature> loadFeatures(BufferedReader reader, Genome genome) {
        return loadFeatures(reader, genome, new GFFCodec(genome));
    }

    public List<htsjdk.tribble.Feature> loadFeatures(BufferedReader reader, Genome genome, GFFCodec codec) {
        String line = null;
        int lineNumber = 0;
        GFFFeatureSource.GFFCombiner combiner = new GFFFeatureSource.GFFCombiner();
        try {
            while ((line = reader.readLine()) != null) {
                lineNumber++;

                if (line.startsWith("#")) {
                    codec.readHeaderLine(line);
                } else {
                    try {
                        Feature f = codec.decode(line);
                        if (f != null) {
                            combiner.addFeature((BasicFeature) f);
                        }
                    } catch (Exception e) {
                        log.error("Error parsing: " + line, e);
                    }
                }
            }


        } catch (IOException ex) {
            log.error("Error reading GFF file", ex);
            if (line != null && lineNumber != 0) {
                throw new ParserException(ex.getMessage(), ex, lineNumber, line);
            } else {
                throw new RuntimeException(ex);
            }
        }

        trackProperties = TrackLoader.getTrackProperties(codec.getHeader());

        //Combine the features
        List<Feature> iFeatures = combiner.combineFeatures();

        FeatureDB.addFeatures(iFeatures, genome);

        return iFeatures;
    }


    public static Set<String> geneParts = new HashSet();

    static {
        geneParts.add("five_prime_UTR");
        geneParts.add("three_prime_UTR");
        geneParts.add("5'-utr");
        geneParts.add("3'-utr");
        geneParts.add("3'-UTR");
        geneParts.add("5'-UTR");
        geneParts.add("5utr");
        geneParts.add("3utr");
        geneParts.add("CDS");
        geneParts.add("cds");
        geneParts.add("exon");
        geneParts.add("coding_exon");
        geneParts.add("intron");
        geneParts.add("transcript");
        geneParts.add("processed_transcript");
        geneParts.add("mrna");
        geneParts.add("mRNA");

    }

    /**
     * Given a GFF File, creates a new GFF file for each type. Any feature type
     * which is part of a "gene" ( {@link SequenceOntology#geneParts} ) are put in the same file,
     * others are put in different files. So features of type "gene", "exon", and "mrna"
     * would go in gene.gff, but features of type "myFeature" would go in myFeature.gff.
     *
     * @param gffFile
     * @param outputDirectory
     * @throws IOException
     */
    public static void splitFileByType(String gffFile, String outputDirectory) throws IOException {

        BufferedReader br = new BufferedReader(new FileReader(gffFile));
        String nextLine;
        String ext = "." + gffFile.substring(gffFile.length() - 4);

        Map<String, PrintWriter> writers = new HashMap();

        while ((nextLine = br.readLine()) != null) {
            nextLine = nextLine.trim();
            if (!nextLine.startsWith("#")) {
                String[] tokens = Globals.tabPattern.split(nextLine.trim().replaceAll("\"", ""), -1);

                String type = tokens[2];
                if (SequenceOntology.geneParts.contains(type)) {
                    type = "gene";
                }
                if (!writers.containsKey(type)) {
                    writers.put(type,
                            new PrintWriter(new FileWriter(new File(outputDirectory, type + ext))));
                }
            }
        }
        br.close();

        br = new BufferedReader(new FileReader(gffFile));
        PrintWriter currentWriter = null;
        while ((nextLine = br.readLine()) != null) {
            nextLine = nextLine.trim();
            if (nextLine.startsWith("#")) {
                for (PrintWriter pw : writers.values()) {
                    pw.println(nextLine);
                }
            } else {
                String[] tokens = Globals.tabPattern.split(nextLine.trim().replaceAll("\"", ""), -1);
                String type = tokens[2];
                if (SequenceOntology.geneParts.contains(type)) {
                    type = "gene";
                }
                currentWriter = writers.get(type);

                if (currentWriter != null) {
                    currentWriter.println(nextLine);
                } else {
                    System.out.println("No writer for: " + type);
                }
            }

        }

        br.close();
        for (PrintWriter pw : writers.values()) {
            pw.close();
        }
    }

    public TrackProperties getTrackProperties() {
        return trackProperties;
    }

    public static void main(String[] args) throws IOException {
        if (args.length < 2) {
            System.out.println("SpitFilesByType <gffFile> <outputDirectory>");
            return;
        }
        splitFileByType(args[0], args[1]);
    }
}
