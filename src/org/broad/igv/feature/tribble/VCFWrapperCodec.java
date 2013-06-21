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
import org.broad.igv.variant.Variant;
import org.broad.igv.variant.vcf.VCFVariant;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.Feature;
import org.broad.tribble.readers.LineReader;
import org.broadinstitute.variant.variantcontext.VariantContext;

/**
 * @author Jim Robinson
 * @date Aug 1, 2011
 */
public class VCFWrapperCodec extends AsciiFeatureCodec<VCFVariant> {

    private static Logger log = Logger.getLogger(Variant.class);

    AsciiFeatureCodec wrappedCodec;
    Genome genome;

    public VCFWrapperCodec(AsciiFeatureCodec wrappedCodec, Genome genome) {
        super(VCFVariant.class);
        this.wrappedCodec = wrappedCodec;
        this.genome = genome;
    }

    @Override
    public Feature decodeLoc(String line) {
        return wrappedCodec.decodeLoc(line);
    }

    @Override
    public VCFVariant decode(String line) {
        // VCFCodec supports completely missing fields (which would simply have a ".")
        // but does not currently support missing only certain elements of a field.
        // IGV is much more permissive.

        VariantContext vc = null;
        try {
            vc = (VariantContext) wrappedCodec.decode(line);
            //The genotype fields are loaded lazily, we force parsing here to
            //catch the exception if necessary
            if (vc != null) vc.getSampleNames();
        } catch (NumberFormatException e) {
            String msg = String.format("NumberFormatException on line: %s \n Attempting to reformat by replacing ,., with ,0,", line);
            log.warn(msg);
            String refLine = line.replaceAll(",\\.", ",0");
            refLine = refLine.replaceAll("\\.,", "0,");
            vc = (VariantContext) wrappedCodec.decode(refLine);
        }


        if (vc == null) {
            return null;
        }
        String chr = genome == null ? vc.getChr() : genome.getChromosomeAlias(vc.getChr());
        return new VCFVariant(vc, chr);

    }

    @Override
    public Object readHeader(LineReader reader) {
        return wrappedCodec.readHeader(reader);
    }

    /**
     * This function returns true iff the File potentialInput can be parsed by this
     * codec.
     * <p/>
     * There is an assumption that there's never a situation where two different Codecs
     * return true for the same file.  If this occurs, the recommendation would be to error out.
     * <p/>
     * Note this function must never throw an error.  All errors should be trapped
     * and false returned.
     *
     * @param path the file to test for parsability with this codec
     * @return true if potentialInput can be parsed, false otherwise
     */
    @Override
    public boolean canDecode(String path) {
        return path.endsWith(".vcf") || path.endsWith(".vcf4") || path.endsWith(".vcf3");
    }
}
