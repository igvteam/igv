/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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

package org.broad.igv.variant.util;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.FeatureCodec;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

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
