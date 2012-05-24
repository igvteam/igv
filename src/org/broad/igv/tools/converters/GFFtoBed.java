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

package org.broad.igv.tools.converters;

import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.GFFParser;
import org.broad.igv.feature.tribble.IGVBEDCodec;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.Feature;

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

        GFFParser parser = new GFFParser(inputFile.getPath());

        BufferedReader reader = null;
        PrintWriter pw = null;
        try {
            reader = ParsingUtils.openBufferedReader(new ResourceLocator(inputFile.getAbsolutePath()));
            List<Feature> features = parser.loadFeatures(reader, null);

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
        File gffFile = new File("test/data/gff/gene.gff3");
        File bedFile = new File("test/data/gff/gene.bed");
        convert(gffFile, bedFile);
    }
}
