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
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.feature;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.exceptions.ParserException;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.renderer.IGVFeatureRenderer;
import org.broad.igv.track.FeatureCollectionSource;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;

import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * @author jrobinso
 */
public abstract class AbstractFeatureParser implements FeatureParser {

    private static Logger log = Logger.getLogger(IGV.class);
    protected int startBase = 0;
    boolean gffTags = false;

    /* An object to collection track properties, if specified in the feature file. */
    protected TrackProperties trackProperties = null;


    public static boolean canParse(String path) {
        return getCodec(path, null) != null;
    }

    /**
     * Return an parser instance appropriate the the file type.  Currently the filename
     * is used to determine file type, this is fragile obviously but what is the alternative?
     */
    public static FeatureParser getInstanceFor(ResourceLocator locator, Genome genome) {

        final String path = locator.getPath();
        return getInstanceFor(path, genome);
    }


    public static FeatureParser getInstanceFor(String path, Genome genome) {
        FeatureCodec codec = getCodec(path, genome);
        if (codec != null) {
            return new FeatureCodecParser(codec, genome);
        } else {
            return null;
        }
    }

    private static FeatureCodec getCodec(String path, Genome genome) {
        String tmp = getStrippedFilename(path);
        return CodecFactory.getCodec(tmp, genome);
    }



    /**
     * Load all features in this file.
     *
     * @param locator
     * @return
     */
    public List<FeatureTrack> loadTracks(ResourceLocator locator, Genome genome) {

        List<Feature> features = loadFeatures(locator, -1);
        if (features.size() == 0) {
            //MessageUtils.showMessage("<html>Warning.  No features were found in " + locator.getPath() + ".<br>Track not loaded.");
        }
        return loadTracks(features, locator, genome, trackProperties);
    }

    private List<FeatureTrack> loadTracks(List<Feature> features, ResourceLocator locator, Genome genome) {
        FeatureCollectionSource source = new FeatureCollectionSource(features, genome);
        FeatureTrack track = new FeatureTrack(locator, source);
        track.setName(locator.getTrackName());
        track.setRendererClass(IGVFeatureRenderer.class);
        track.setHeight(45);

        //Nothing writes to trackType as far as I can tell -Jacob S
        //if (trackType != null) {
        //    track.setTrackType(trackType);
        //}
        if (trackProperties != null) {
            track.setProperties(trackProperties);
        }

        List<FeatureTrack> tracks = new ArrayList();
        tracks.add(track);
        return tracks;
    }

    public static List<FeatureTrack> loadTracks(Iterable<Feature> features, ResourceLocator locator, Genome genome,
                                                TrackProperties trackProperties) {


        FeatureCollectionSource source = new FeatureCollectionSource(features, genome);

        //Load into FeatureDB for searching
        if (IGV.hasInstance() || Globals.isTesting()) {
            Set<String> chrs = source.getChrs();
            for (String chr : chrs) {
                List<Feature> feats = source.getFeatures(chr);
                for (Feature f : feats) {
                    if (f instanceof NamedFeature) {
                        FeatureDB.addFeature((NamedFeature) f);
                    }
                }
            }
        }

        FeatureTrack track = new FeatureTrack(locator, source);
        track.setName(locator.getTrackName());
        track.setRendererClass(IGVFeatureRenderer.class);
        track.setHeight(45);

        if (trackProperties != null) {
            track.setProperties(trackProperties);
        }

        List<FeatureTrack> tracks = new ArrayList();
        tracks.add(track);

        return tracks;
    }

    /**
     * Parse a limited number of lines in this file and return a list of features found.
     *
     * @param locator
     * @param maxLines
     * @return
     */
    public List<org.broad.tribble.Feature> loadFeatures(ResourceLocator locator, int maxLines) {

        BufferedReader reader = null;
        try {
            reader = ParsingUtils.openBufferedReader(locator);
            return loadFeatures(reader, maxLines);
        } catch (IOException e) {
            throw new RuntimeException(e);
        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e) {

                }
            }
        }
    }

    /**
     * Method description
     *
     * @param reader
     * @return
     */
    public List<org.broad.tribble.Feature> loadFeatures(BufferedReader reader, Genome genome) {
        return loadFeatures(reader, -1);
    }

    /**
     * Load all features in this file.
     *
     * @param reader
     * @param maxLines
     * @return
     */
    public List<org.broad.tribble.Feature> loadFeatures(BufferedReader reader, int maxLines) {

        List<org.broad.tribble.Feature> features = new ArrayList<org.broad.tribble.Feature>();
        String nextLine = null;

        int nLines = 0;
        try {
            while ((nextLine = reader.readLine()) != null) {
                nextLine = nextLine.trim();
                if (nextLine.length() == 0) continue;
                nLines++;
                if ((maxLines > 0) && (nLines > maxLines)) {
                    break;
                }

                try {
                    if (nextLine.startsWith("#")) {
                        if (nextLine.startsWith("#type")) {
                            String[] tokens = Globals.equalPattern.split(nextLine);
                            if (tokens.length > 1) {
                                try {
                                    // TODO: type is not currently used, is there any reason to keep this?
                                    TrackType type = TrackType.valueOf(tokens[1]);
                                } catch (Exception e) {
                                    log.error("Error converting track type: " + tokens[1]);
                                }
                            }
                        } else if (nextLine.startsWith("#track")) {
                            TrackProperties tp = new TrackProperties();
                            ParsingUtils.parseTrackLine(nextLine, tp);
                            setTrackProperties(tp);
                            if (tp.isGffTags()) {
                                gffTags = true;
                            }
                        } else if (nextLine.startsWith("#coords")) {
                            try {
                                String[] tokens = Globals.equalPattern.split(nextLine);
                                startBase = Integer.parseInt(tokens[1]);
                            } catch (Exception e) {
                                log.error("Error parsing coords line: " + nextLine, e);
                            }

                        } else if (nextLine.startsWith("#gffTags")) {
                            gffTags = true;
                        }
                    } else {
                        Feature feature = parseLine(nextLine);
                        if (feature != null) {
                            features.add(feature);
                        }
                    }

                } catch (NumberFormatException e) {

                    // Expected condition -- for example comments.  don't log as it slows down
                    // the parsing and is not useful information.
                }
            }
        } catch (java.io.EOFException e) {

            // This exception is due to a known bug with java zip library.  Not
            // in general a real error, and nothing we can do about it in any
            // event.
            return features;
        } catch (Exception e) {
            if (nextLine != null && nLines != 0) {
                throw new ParserException(e.getMessage(), e, nLines, nextLine);
            } else {
                throw new RuntimeException(e);
            }
        }

        // TODO -- why is this test here?  This will break igvtools processing of expression files
        //if (IGV.hasInstance() || Globals.isTesting()) {
            FeatureDB.addFeatures(features);
        //}
        return features;
    }

    abstract protected Feature parseLine(String nextLine);

    protected static String getStrippedFilename(String filename) {
        String tmp = filename.toLowerCase();

//      String off common, non-informative extensions .gz, .tab, .txt, and .csv
        if (filename.endsWith(".gz")) {
            tmp = tmp.substring(0, tmp.length() - 3);
        }
        if (tmp.endsWith(".tab") || tmp.endsWith(".txt") || tmp.endsWith(".csv")) {
            tmp = tmp.substring(0, tmp.length() - 4);
        }
        return tmp;
    }

    /**
     * Convenience method.  Write a list of features out as a BED file
     *
     * @param features
     * @param outputfile
     */
    public static void dumpFeatures(List<IGVFeature> features, String outputfile) {

        PrintWriter pw = null;

        try {
            pw = new PrintWriter(new FileWriter(outputfile));
            pw.println("Header row");
            for (IGVFeature gene : features) {
                pw.print(gene.getName() + "\t");
                pw.print(gene.getIdentifier() + "\t");
                pw.print(gene.getChr() + "\t");
                if (gene.getStrand() == Strand.POSITIVE) {
                    pw.print("+\t");
                } else if (gene.getStrand() == Strand.NEGATIVE) {
                    pw.print("-\t");
                } else {
                    pw.print(" \t");
                }
                pw.print(gene.getStart() + "\t");
                pw.print(gene.getEnd() + "\t");

                List<Exon> regions = gene.getExons();
                pw.print(regions.size() + "\t");
                for (Exon exon : regions) {
                    pw.print(exon.getStart() + ",");
                }
                pw.print("\t");
                for (Exon exon : regions) {
                    pw.print(exon.getEnd() + ",");
                }
                pw.println();

            }
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (pw != null) {
                pw.close();
            }
        }
    }

    protected void setTrackProperties(TrackProperties trackProperties) {
        this.trackProperties = trackProperties;
    }

    public TrackProperties getTrackProperties() {
        return trackProperties;
    }

}
