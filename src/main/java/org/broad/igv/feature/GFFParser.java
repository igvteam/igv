package org.broad.igv.feature;

import htsjdk.tribble.Feature;
import org.broad.igv.logging.*;
import org.broad.igv.exceptions.ParserException;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.gff.GFFCombiner;
import org.broad.igv.feature.gff.GFFFeatureSource;
import org.broad.igv.feature.tribble.GFFCodec;
import org.broad.igv.track.*;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @deprecated use org.broad.igv.feature.gff.GFFFeatureSource
 * User: jrobinso
 */

@Deprecated
public class GFFParser implements FeatureParser {

    static Logger log = LogManager.getLogger(GFFParser.class);

    private TrackProperties trackProperties = null;

    public List<htsjdk.tribble.Feature> loadFeatures(BufferedReader reader, Genome genome) {
        return loadFeatures(reader, genome, new GFFCodec(genome, null));
    }

    public List<htsjdk.tribble.Feature> loadFeatures(BufferedReader reader, Genome genome, GFFCodec codec) {
        String line = null;
        int lineNumber = 0;
        GFFCombiner combiner = GFFFeatureSource.getCombiner(codec.getVersion());
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

        if(genome != null) {
            genome.getFeatureDB().addFeatures(iFeatures);
        }

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
