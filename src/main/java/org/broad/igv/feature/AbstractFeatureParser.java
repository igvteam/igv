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
import org.broad.igv.feature.tribble.FeatureFileHeader;
import org.broad.igv.track.FeatureCollectionSource;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;

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


    public static boolean canParse(ResourceLocator locator) {
        return CodecFactory.getCodec(locator, null) != null;
    }

    /**
     * Return an parser instance appropriate the the file type.  Currently the filename
     * is used to determine file type, this is fragile obviously but what is the alternative?
     */
    public static FeatureParser getInstanceFor(ResourceLocator locator, Genome genome) {
        FeatureCodec codec = CodecFactory.getCodec(locator, genome);
        if (codec != null && codec instanceof AsciiFeatureCodec) {
            return new FeatureCodecParser((AsciiFeatureCodec) codec, genome);
        } else {
            return null;
        }
    }

    /**
     *
     * @param reader
     * @return
     */
    public List<htsjdk.tribble.Feature> loadFeatures(BufferedReader reader, Genome genome) {
        return loadFeatures(reader, genome, -1);
    }

    /**
     * Load all features in this file.
     *
     * @param reader
     * @param genome
     * @param maxLines
     * @return
     */
    public List<htsjdk.tribble.Feature> loadFeatures(BufferedReader reader, Genome genome, int maxLines) {

        List<htsjdk.tribble.Feature> features = new ArrayList<htsjdk.tribble.Feature>();
        String nextLine = null;

        int maxLogErrors = 10;
        int nErrors = 0;
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
                        } else if (nextLine.startsWith("#track") || nextLine.startsWith("track")) {
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
                    if (nErrors < maxLogErrors) {
                        log.error("Number format error parsing line: " + nextLine, e);
                    }
                    nErrors++;
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
        FeatureDB.addFeatures(features, genome);
        //}
        return features;
    }

    abstract protected Feature parseLine(String nextLine);

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
