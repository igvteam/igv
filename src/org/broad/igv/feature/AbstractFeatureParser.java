/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.feature;

import org.apache.log4j.Logger;
import org.broad.igv.exceptions.ParserException;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.renderer.IGVFeatureRenderer;
import org.broad.igv.track.FeatureCollectionSource;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.IGVMainFrame;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * @author jrobinso
 */
public abstract class AbstractFeatureParser implements FeatureParser {

    private static Logger log = Logger.getLogger(IGVMainFrame.class);
    private String filePath;
    protected Genome genome;
    protected int startBase = 0;

    /* An object to collection track properties, if specified in the feature file. */
    protected TrackProperties trackProperties = null;

    TrackType trackType;

    public AbstractFeatureParser() {
        genome = GenomeManager.getInstance().getCurrentGenome();
    }

    /**
     * Return an parser instance appropriate the the file type.  Currently the filename
     * is used to determine file type, this is fragile obviously but what is the alternative?
     */
    public static FeatureParser getInstanceFor(ResourceLocator locator) {

        final String path = locator.getPath();
        return getInstanceFor(path);
    }

    public static FeatureParser getInstanceFor(String path) {
        String tmp = getStrippedFilename(path);

        if (tmp.endsWith("bed") || tmp.endsWith("map")) {
            return new BEDFileParser();
        } else if (tmp.contains("refflat")) {
            return new UCSCGeneTableParser(UCSCGeneTableParser.Type.REFFLAT);
        } else if (tmp.contains("genepred") || tmp.contains("ensgene") || tmp.contains("refgene")) {
            return new UCSCGeneTableParser(UCSCGeneTableParser.Type.GENEPRED);
        } else if (tmp.contains("ucscgene")) {
            return new UCSCGeneTableParser(UCSCGeneTableParser.Type.UCSCGENE);
        } else if (tmp.endsWith("gtf") || tmp.endsWith("gff") || tmp.endsWith("gff3")) {
            return new GFFParser(path);
        } else if (tmp.endsWith("embl")) {
            return new EmblFeatureTableParser();
        } else {
            return null;
        }
    }


    /**
     * Return true if the file represented by the locator contains feature.  This method
     * returns true by default.  It can be overriden by subclasses representing ambiguous
     * file content.
     *
     * @param locator
     * @return true if a feature file
     */
    public boolean isFeatureFile(ResourceLocator locator) {
        return true;

    }

    /**
     * Load all features in this file.
     *
     * @param locator
     * @return
     */
    public List<FeatureTrack> loadTracks(ResourceLocator locator) {

        List<org.broad.tribble.Feature> features = loadFeatures(locator, -1);
        if(features.size() == 0) {
            throw new DataLoadException("No features found in this file with locations mapped to the current genome", locator.getPath());
        }
        FeatureCollectionSource source = new FeatureCollectionSource(features);
        FeatureTrack track = new FeatureTrack(locator, source);
        track.setName(locator.getTrackName());
        track.setRendererClass(IGVFeatureRenderer.class);
        track.setMinimumHeight(35);
        track.setHeight(45);
        //track.setRendererClass(GeneTrackRenderer.class);
        if (trackType != null) {
            track.setTrackType(trackType);
        }
        if (trackProperties != null) {
            track.setTrackProperties(trackProperties);
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
    protected List<org.broad.tribble.Feature> loadFeatures(ResourceLocator locator, int maxLines) {

        // File file = new File(locator.getPath());
        int nLines = 0;
        AsciiLineReader reader = null;
        try {
            reader = ParsingUtils.openAsciiReader(locator);
            return loadFeatures(reader, maxLines);
        } catch (IOException e) {
            throw new RuntimeException(e);
        } finally {
            if (reader != null) {
                reader.close();
            }
        }
    }

    /**
     * Method description
     *
     * @param reader
     * @return
     */
    public List<org.broad.tribble.Feature> loadFeatures(AsciiLineReader reader) {
        return loadFeatures(reader, -1);
    }

    /**
     * Load all features in this file.
     *
     * @param reader
     * @param maxLines
     * @return
     */
    public List<org.broad.tribble.Feature> loadFeatures(AsciiLineReader reader, int maxLines) {

        List<org.broad.tribble.Feature> features = new ArrayList<org.broad.tribble.Feature>();
        String nextLine = null;

        try {
            int nLines = 0;
            while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
                nextLine = nextLine.trim();
                nLines++;
                if ((maxLines > 0) && (nLines > maxLines)) {
                    break;
                }

                try {
                    if (nextLine.startsWith("#")) {
                        if (nextLine.startsWith("#type")) {
                            String[] tokens = nextLine.split("=");
                            if (tokens.length > 1) {
                                try {
                                    TrackType type = TrackType.valueOf(tokens[1]);
                                } catch (Exception e) {
                                    log.error("Error converting track type: " + tokens[1]);
                                }
                            }
                        } else if (nextLine.startsWith("#track")) {
                            TrackProperties tp = new TrackProperties();
                            ParsingUtils.parseTrackLine(nextLine, tp);
                            setTrackProperties(tp);
                        } else if(nextLine.startsWith("#coords")) {
                            try {
                                String[] tokens = nextLine.split("=");
                                startBase = Integer.parseInt(tokens[1]);
                            }
                            catch(Exception e) {
                                log.error("Error parsing coords line: " + nextLine, e);
                            }
                            
                        }
                    } else {
                        IGVFeature feature = parseLine(nextLine);
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
            if (nextLine != null && reader.getCurrentLineNumber() != 0) {
                throw new ParserException(e.getMessage(), e, reader.getCurrentLineNumber(), nextLine);
            } else {
                throw new RuntimeException(e);
            }
        }

        //parsingComplete(features);
        FeatureDB.addFeatures(features);
        return features;
    }

    protected void parsingComplete(List<IGVFeature> features) {

        // Default -- do nothing.
    }

    abstract protected IGVFeature parseLine(String nextLine);

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
