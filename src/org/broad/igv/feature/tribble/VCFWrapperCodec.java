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

import org.apache.log4j.Logger;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.variant.Variant;
import org.broad.igv.variant.vcf.VCFVariant;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.Feature;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;

import java.io.IOException;

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
    public Feature decodeLoc(LineIterator iterator) throws IOException{
        return wrappedCodec.decodeLoc(iterator);
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
    public Object readActualHeader(LineIterator reader) {
        return wrappedCodec.readActualHeader(reader);
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
