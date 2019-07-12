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

package org.broad.igv.tools.converters;

import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.GFFParser;
import org.broad.igv.feature.tribble.GFFCodec;
import org.broad.igv.feature.tribble.IGVBEDCodec;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
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
            GFFCodec gffCodec = new GFFCodec(version, null);
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
