package org.igv.tools.converters;

import org.igv.feature.BasicFeature;
import org.igv.feature.GFFParser;
import org.igv.feature.tribble.GFFCodec;
import org.igv.feature.tribble.IGVBEDCodec;
import org.igv.util.ParsingUtils;
import org.igv.util.ResourceLocator;
import htsjdk.tribble.Feature;

import java.io.*;
import java.util.List;

/**
 * Convert a GFF to a BED file, keeping the column 9 attribute tags.
 *
 * @author Jim Robinson
 * @date 2/5/12
 */
public class GFFtoBed {


    public static void convert(File inputFile, File outputFile) {

        GFFParser parser = new GFFParser();

        BufferedReader reader = null;
        PrintWriter pw = null;
        try {
            GFFCodec.Version version = inputFile.getPath().endsWith(".gff3") ? GFFCodec.Version.GFF3 : GFFCodec.Version.GFF2;
            GFFCodec gffCodec = new GFFCodec(version, null, null);
            reader = ParsingUtils.openBufferedReader(new ResourceLocator(inputFile.getAbsolutePath()));
            List<Feature> features = parser.loadFeatures(reader, null, gffCodec);

            IGVBEDCodec codec = new IGVBEDCodec();
            codec.setGffTags(true);

            pw = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));
            pw.println("#gffTags");
            for(Feature feature : features) {

                BasicFeature bf = (BasicFeature) feature;
                String line = codec.encode(bf);
                pw.println(line);
            }


        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e) {

                }
            }
            try {
                if (pw != null) pw.close();
            } catch (Exception e) {

            }
        }


    }


    public static void main(String [] args) {
        File gffFile = new File("/Users/jrobinso/projects/igv-js/igv/assets/hg19/gencode.v18.transcripts.patched_contigs.gtf");
        File bedFile = new File("/Users/jrobinso/projects/igv-js/igv/assets/hg19/gencode.v18.transcripts.patched_contigs.bed");
        convert(gffFile, bedFile);
    }
}
