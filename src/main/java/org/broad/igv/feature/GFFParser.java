/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2018 Broad Institute
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

import htsjdk.tribble.Feature;
import org.apache.log4j.Logger;
import org.broad.igv.exceptions.ParserException;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.GFFCodec;
import org.broad.igv.track.*;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @deprecated use org.broad.igv.track.GFFFeatureSource
 * User: jrobinso
 */

@Deprecated
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

    public TrackProperties getTrackProperties() {
        return trackProperties;
    }

    public static void main(String[] args) throws IOException {
        if (args.length < 2) {
            System.out.println("SpitFilesByType <gffFile> <outputDirectory>");
            return;
        }
        FeatureFileUtils.splitGffFileByType(args[0], args[1]);
    }
}
