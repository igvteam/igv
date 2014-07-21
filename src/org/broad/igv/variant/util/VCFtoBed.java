/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.variant.util;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.CodecFactory;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureCodec;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;

/**
 * Utility class to convert a VCF to a BED file.  Genotypes if any are lost in this converstion
 *
 * @author jrobinso
 * @date Jan 22, 2011
 */
public class VCFtoBed {


    public static void main(String[] args) throws IOException {
        convert(args[0], args[1]);
    }

    public static void convert(String vcfFile, String bedFile) throws IOException {

        AbstractFeatureReader basicReader = null;
        PrintWriter writer = null;
        try {
            writer = new PrintWriter(new BufferedWriter(new FileWriter(bedFile)));
            Genome genome = null; // Don't do alias converting
            FeatureCodec codec = CodecFactory.getCodec(vcfFile, genome);
            boolean isVCF = codec.getClass().isAssignableFrom(VCFCodec.class);
            basicReader = AbstractFeatureReader.getFeatureReader(vcfFile, codec, true);

            Iterator<VariantContext> iter = basicReader.iterator();

            while (iter.hasNext()) {

                VariantContext vc = iter.next();
                String chr = vc.getChr();
                if (!chr.startsWith("chr")) {
                    chr = "chr" + chr;
                }

                int start = vc.getStart() - 1;
                int end = vc.getEnd();
                String id = vc.getID();
                if(id == null) {
                    id = ".";
                }
                String af = vc.getAttributeAsString("AF", "");

                writer.println(chr + "\t" + start + "\t" + end + "\t" + id + "\t" + af);

            }
        }
        finally {
            if (basicReader != null)
                basicReader.close();
            if (writer != null)
                writer.close();
        }


    }


}
