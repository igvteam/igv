/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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
import org.broad.tribble.Feature;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * User: jrobinso
 */


public class GFFParser implements FeatureParser {

    static Logger log = Logger.getLogger(GFFParser.class);

    private TrackProperties trackProperties = null;

    public GFFParser(String path) {
    }

    public List<FeatureTrack> loadTracks(ResourceLocator locator, Genome genome) {

        BufferedReader reader = null;
        try {
            reader = ParsingUtils.openBufferedReader(locator);
            GFFCodec.Version version = locator.getPath().endsWith(".gff3") ? GFFCodec.Version.GFF3 : GFFCodec.Version.GFF2;
            List<org.broad.tribble.Feature> features = loadFeatures(reader, genome);

            FeatureTrack track = new FeatureTrack(locator, new FeatureCollectionSource(features, genome));
            track.setName(locator.getTrackName());
            track.setRendererClass(IGVFeatureRenderer.class);
            track.setMinimumHeight(35);
            track.setHeight(45);
            track.setRendererClass(GeneTrackRenderer.class);

            if (trackProperties != null) {
                track.setProperties(trackProperties);
                track.setName(trackProperties.getName());
            }

            List<FeatureTrack> tracks = new ArrayList();
            tracks.add(track);
            return tracks;

        } catch (IOException ex) {
            log.error(ex);
            throw new RuntimeException(ex);

        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e) {

                }

            }
        }
    }

    public List<org.broad.tribble.Feature> loadFeatures(BufferedReader reader, Genome genome) {
        List<org.broad.tribble.Feature> features = new ArrayList();
        String line = null;
        int lineNumber = 0;
        GFFCodec codec = new GFFCodec(genome);
        try {


            while ((line = reader.readLine()) != null) {
                lineNumber++;

                if (line.startsWith("#")) {
                    codec.readHeaderLine(line);
                } else {
                    try {
                        Feature f = codec.decode(line);
                        if (f != null) {
                            features.add(f);
                        }
                    } catch (Exception e) {
                        e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
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
        List<Feature> iFeatures = (new GFFFeatureSource.GFFCombiner()).combineFeatures(features.iterator());

        if (IGV.hasInstance()) {
            FeatureDB.addFeatures(iFeatures);
        }

        return iFeatures;
    }


    /**
     * Given a GFF File, creates a new GFF file for each type. Any feature type
     * which is part of a "gene" ( {@link GFFCodec#geneParts} ) are put in the same file,
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
                if (GFFCodec.geneParts.contains(type)) {
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
                if (GFFCodec.geneParts.contains(type)) {
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
