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

package org.broad.igv.feature.tribble;

import org.apache.log4j.Logger;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.peaks.PeakCodec;
import org.broad.igv.util.ParsingUtils;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.util.BlockCompressedInputStream;
import org.broadinstitute.sting.utils.codecs.vcf.VCF3Codec;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * A factory class for Tribble codecs.  implements a single, static, public method to return a codec given a
 * path to a feature file  (bed, gff, vcf, etc).
 */
public class CodecFactory {

    private static Logger log = Logger.getLogger(CodecFactory.class);

    public static final List<String> validExtensions = new ArrayList<String>(15);

    static {
        validExtensions.addAll(Arrays.asList("vcf4", "vcf", "bed", "refflat", "genepred", "ensgene", "refgene", "ucscgene",
                "repmask", "gff3", "gvf", "gff", "gtf", "psl"));
    }

    /**
     * Return a tribble codec to decode the supplied file, or null if not found.
     *
     * @param path the path (file or URL) to the feature file
     */
    public static AsciiFeatureCodec getCodec(String path, Genome genome) {

        String fn = path.toLowerCase();
        if (fn.endsWith(".gz")) {
            int l = fn.length() - 3;
            fn = fn.substring(0, l);
        }
        if (fn.endsWith(".txt")) {
            int l = fn.length() - 4;
            fn = fn.substring(0, l);
        }

        if (fn.endsWith(".vcf3")) {
            return new VCFWrapperCodec(new VCF3Codec(), genome);
        }
        if (fn.endsWith(".vcf4")) {
            return new VCFWrapperCodec(new VCFCodec(), genome);
        } else if (fn.endsWith(".vcf")) {
            return new VCFWrapperCodec(getVCFCodec(path), genome);
        } else if (fn.endsWith(".bed")) {
            return new IGVBEDCodec(genome);
        } else if (fn.contains("refflat")) {
            return new UCSCGeneTableCodec(UCSCGeneTableCodec.Type.REFFLAT, genome);
        } else if (fn.contains("genepred") || fn.contains("ensgene") || fn.contains("refgene")) {
            return new UCSCGeneTableCodec(UCSCGeneTableCodec.Type.GENEPRED, genome);
        } else if (fn.contains("ucscgene")) {
            return new UCSCGeneTableCodec(UCSCGeneTableCodec.Type.UCSCGENE, genome);
        } else if (fn.endsWith(".rmask") || (fn.endsWith(".repmask"))) {
            return new REPMaskCodec(genome);
        } else if (fn.endsWith(".gff3") || fn.endsWith(".gvf")) {
            return new GFFCodec(GFFCodec.Version.GFF3, genome);
        } else if (fn.endsWith(".gff") || fn.endsWith(".gtf")) {
            return new GFFCodec(genome);
            //} else if (fn.endsWith(".sam")) {
            //return new SAMCodec();
        } else if (fn.endsWith(".psl") || fn.endsWith(".pslx")) {
            return new PSLCodec(genome);

        }
//   Mutation codec disabled until we deal with the multiple sample problem
//        else if (fn.endsWith(".mut") || (fn.endsWith(".maf") && MUTCodec.isMutationAnnotationFile(path))) {
//            return new MUTCodec(path, genome);
//        }
        else if (fn.endsWith(".narrowpeak") || fn.endsWith(".broadpeak")) {
            return new EncodePeakCodec(genome);
        } else if (fn.endsWith(".peak")) {
            return new PeakCodec(genome);

        } else {
            return null;
        }

    }


    /**
     * Return the appropriate VCFCodec based on the version tag.
     * <p/>
     * e.g.  ##fileformat=VCFv4.1
     *
     * @param path Path to the VCF file.  Can be a file path, or URL
     * @return
     */
    private static AsciiFeatureCodec getVCFCodec(String path) {

        BufferedReader reader = null;

        try {
            // If the file ends with ".gz" assume it is a tabix indexed file
            if (path.toLowerCase().endsWith(".gz")) {
                reader = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(
                        org.broad.tribble.util.ParsingUtils.openInputStream(path))));
            } else {
                reader = ParsingUtils.openBufferedReader(path);
            }
            // Look for fileformat directive.  This should be the first line, but just in case check the first 20
            int lineCount = 0;
            String formatLine;
            while ((formatLine = reader.readLine()) != null && lineCount < 20) {
                if (formatLine.toLowerCase().startsWith("##fileformat")) {
                    String[] tmp = formatLine.split("=");
                    if (tmp.length > 1) {
                        String version = tmp[1].toLowerCase();
                        if (version.startsWith("vcfv3")) {
                            return new VCF3Codec();
                        } else {
                            return new VCFCodec();
                        }
                    }
                }
                lineCount++;
            }

        } catch (IOException e) {
            log.error("Error checking VCF Version");

        } finally {
            if (reader != null) try {
                reader.close();
            } catch (IOException e) {

            }
        }
        // Should never get here, but as a last resort assume this is a VCF 4.x file.
        return new VCFCodec();
    }
}
