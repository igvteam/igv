package org.igv.feature.tribble;

import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.FeatureCodec;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.vcf.VCF3Codec;
import htsjdk.variant.vcf.VCFCodec;
import org.igv.data.cufflinks.FPKMTrackingCodec;
import org.igv.feature.FeatureType;
import org.igv.feature.dsi.DSICodec;
import org.igv.feature.genome.Genome;
import org.igv.bedpe.InteractCodec;
import org.igv.gwas.EQTLCodec;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.util.FileUtils;
import org.igv.util.HttpUtils;
import org.igv.util.ParsingUtils;
import org.igv.util.ResourceLocator;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;

/**
 * A factory class for Tribble codecs.  implements a single, static, public method to return a codec given a
 * path to a feature file  (bed, gff, vcf, etc).
 */
public class CodecFactory {

    private static Logger log = LogManager.getLogger(CodecFactory.class);

    public static final List<String> validExtensions = new ArrayList<String>(15);

    public static String ucscSNP = "snp[0-9][0-9][0-9]";

    /**
     * @param path
     * @param genome This won't handle URLs with query strings properly for all codecs
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
        String format = locator.getFormat();

        switch (format) {
            case "vcf3":
                return new VCFWrapperCodec(new VCF3Codec(), genome);
            case "vcf4":
                return new VCFWrapperCodec(new VCFCodec(), genome);
            case "vcf":
            case "gvcf":
                return new VCFWrapperCodec(getVCFCodec(locator), genome);
            case "bcf":
                return new BCF2WrapperCodec(new BCF2Codec(), genome);
            case "bed":
                return new IGVBEDCodec(genome);
            case "junctions":
                final IGVBEDCodec codec = new IGVBEDCodec(genome);
                codec.setFeatureType(FeatureType.SPLICE_JUNCTION);
                return codec;
            case "gappedpeak":
                return new IGVBEDCodec(genome, FeatureType.GAPPED_PEAK);
            case "dgv":
                return new DGVCodec(genome);
            case "rmask":
            case "repmask":
                return new REPMaskCodec(genome);
            case "gff3":
            case "gvf":
                return new GFFCodec(GFFCodec.Version.GFF3, genome, locator.getFilterTypes());
            case "gff":
                return new GFFCodec(genome, locator.getFilterTypes());
            case "gtf":
                return new GFFCodec(GFFCodec.Version.GTF, genome, locator.getFilterTypes());
            case "psl":
            case "pslx":
                return new PSLCodec(genome);
            case "narrowpeak":
            case "broadpeak":
            case "regionpeak":
                return new EncodePeakCodec(genome);
            case "snp":
            case "ucscsnp":
                return new UCSCSnpCodec(genome);
            case "eqtl":
                return new EQTLCodec(genome);
            case "fpkm_tracking":
                return new FPKMTrackingCodec(path);
            case "dsi":
                return new DSICodec(genome);
            case "paf":
                return new PAFCodec(path, genome);
            case "interval_list":
                return new IntervalListCodec(genome);
            case "refflat":
                return new UCSCGeneTableCodec(UCSCGeneTableCodec.Type.REFFLAT, genome);
            case "refgene":
                return new UCSCGeneTableCodec(UCSCGeneTableCodec.Type.GENEPRED, genome);
            case "ucscgene":
                return new UCSCGeneTableCodec(UCSCGeneTableCodec.Type.UCSCGENE, genome);
            case "genepredext":
                return new UCSCGeneTableCodec(UCSCGeneTableCodec.Type.GENEPRED_EXT, genome);
            case "bedmethyl":
                return new IGVBEDCodec(genome, FeatureType.BED_METHYL);
//            case "interact":
//                return new InteractCodec(genome, FeatureType.INTERACT);

            default:
                if (MUTCodec.isMutationAnnotationFile(locator)) {
                    return new MUTCodec(path, genome);
                } else {
                    return null;
                }
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

        if (FileUtils.isRemote(path)) {
            try {
                path = HttpUtils.createURL(path).toString();
            } catch (MalformedURLException e) {
                log.error("Eror translating url", e);
            }
        }

        BufferedReader reader = null;

        try {
            if (locator.getURLPath().toLowerCase().endsWith(".gz")) {
                reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(
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

        String format = locator.getFormat();

        // The vcf extension is for performance, instantiating a vcf codec is expensive
        return format.equals("vcf") || format.equals("bcf") || getCodec(locator, genome) != null;


    }
}
