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

package org.broad.igv.feature.tribble;

import htsjdk.samtools.util.BlockCompressedInputStream;
import org.apache.log4j.Logger;
import org.broad.igv.data.cufflinks.FPKMTrackingCodec;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.gwas.EQTLCodec;
import org.broad.igv.peaks.PeakCodec;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.FeatureCodec;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.vcf.VCF3Codec;
import htsjdk.variant.vcf.VCFCodec;

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
        validExtensions.addAll(Arrays.asList("vcf4", "vcf", "bed", "refflat", "genepred", "ensgene", "refgene", "ucscgene", "repmask", "gff3", "gvf", "gff", "gtf", "psl", "mut", "maf"));
    }

    /**
     * @param path
     * @param genome
     * @return
     * @deprecated Use {@link #getCodec(org.broad.igv.util.ResourceLocator, org.broad.igv.feature.genome.Genome)}
     * This won't handle URLs with query strings properly for all codecs
     */
    public static FeatureCodec getCodec(String path, Genome genome) {
        return getCodec(new ResourceLocator(path), genome);
    }

    /**
     * Return a tribble codec to decode the supplied file, or null if not found.
     *
     * @param locator the ResourceLocator (file or URL) to the feature file
     */
    public static FeatureCodec getCodec(ResourceLocator locator, Genome genome) {

        String path = locator.getPath();
        String fn = locator.getTypeString().toLowerCase();

        if (fn.endsWith(".vcf3")) {
            return new VCFWrapperCodec(new VCF3Codec(), genome);
        }
        if (fn.endsWith(".vcf4")) {
            return new VCFWrapperCodec(new VCFCodec(), genome);
        } else if (fn.endsWith(".vcf")) {
            return new VCFWrapperCodec(getVCFCodec(locator), genome);
        } else if (fn.endsWith(".bcf")) {
            return new BCF2WrapperCodec(new BCF2Codec(), genome);
        } else if (fn.endsWith(".bed")) {
            final IGVBEDCodec codec = new IGVBEDCodec(genome);
            if (fn.endsWith("junctions.bed")) {
                codec.setSpliceJunctions(true);
            }
            return codec;
        } else if (fn.endsWith(".gappedpeak")) {
            return new IGVBEDCodec(genome, IGVBEDCodec.FeatureType.GAPPED_PEAK);
        }else if (fn.endsWith(".dgv")) {
            return new DGVCodec(genome);
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
        } else if (MUTCodec.isMutationAnnotationFile(locator)) {
            return new MUTCodec(path, genome);
        } else if (fn.endsWith(".narrowpeak") || fn.endsWith(".broadpeak") ) {
            return new EncodePeakCodec(genome);
        } else if (fn.endsWith(".peak")) {
            return new PeakCodec(genome);
        } else if (fn.endsWith(".snp")) {
            return new UCSCSnpCodec(genome);
        } else if (fn.endsWith(".eqtl")) {
            return new EQTLCodec(genome);
        } else if (fn.endsWith("fpkm_tracking")) {
            return new FPKMTrackingCodec(path);
            //} else if (fn.endsWith("gene_exp.diff") || fn.endsWith("cds_exp.diff")) {
            //    return new ExpDiffCodec(path);
        } else {
            return null;
        }

    }


    /**
     * Return the appropriate VCFCodec based on the version tag.
     * <p/>
     * e.g.  ##fileformat=VCFv4.1
     *
     * @param locator
     * @return
     */
    private static AsciiFeatureCodec getVCFCodec(ResourceLocator locator) {

        String path = locator.getPath();

        BufferedReader reader = null;

        try {
            // If the file ends with ".gz" assume it is a tabix indexed file
            if (locator.getURLPath().toLowerCase().endsWith(".gz")) {
                // NOTE:  MUST USE THE PICARD VERSION OF ParsingUtils.  The IGV version will return a gzip stream.
                reader = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(
                        htsjdk.tribble.util.ParsingUtils.openInputStream(path))));
            } else {
                reader = ParsingUtils.openBufferedReader(path);
            }
            // Look for fileformat directive.  This should be the first line, but just in case check the first 20
            int lineCount = 0;
            String formatLine;
            while ((formatLine = reader.readLine()) != null && lineCount < 20) {
                if (formatLine.toLowerCase().startsWith("##fileformat") ||
                        formatLine.toLowerCase().startsWith("##format")) {
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

    /**
     * Return true if a file represented by "path" is indexable.  This method is an optimization, we could just look
     * for the index but that is expensive to do for remote resources.  All tribble indexable extensions should be
     * listed here.
     *
     * @param locator
     * @param genome
     * @return
     */
    public static boolean hasCodec(ResourceLocator locator, Genome genome) {

        String fn = locator.getTypeString();
        if (fn.endsWith(".gz")) {
            int l = fn.length() - 3;
            fn = fn.substring(0, l);
        }
        // The vcf extension is for performance, it doesn't matter which codec is returned all vcf files
        // are indexable.
        return fn.endsWith(".vcf") || fn.endsWith(".bcf") || getCodec(locator, genome) != null;


    }
}
