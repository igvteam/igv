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

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.variant.vcf.VCFVariant;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.FeatureCodecHeader;
import org.broad.tribble.readers.PositionalBufferedStream;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.io.IOException;

/**
 * @author Jacob Silterra
 * @date 2013-Jun-14
 */
public class BCF2WrapperCodec implements FeatureCodec<VCFVariant> {

    FeatureCodec wrappedCodec;
    Genome genome;

    public BCF2WrapperCodec(FeatureCodec wrappedCodec, Genome genome) {
        this.wrappedCodec = wrappedCodec;
        this.genome = genome;
    }

    @Override
    public VCFVariant decode(PositionalBufferedStream stream) throws IOException {
        VariantContext vc = (VariantContext) wrappedCodec.decode(stream);
        if (vc == null) {
            return null;
        }
        String chr = genome == null ? vc.getChr() : genome.getChromosomeAlias(vc.getChr());
        return new VCFVariant(vc, chr);

    }

    @Override
    public Feature decodeLoc(PositionalBufferedStream stream) throws IOException {
        return this.wrappedCodec.decodeLoc(stream);
    }

    @Override
    public FeatureCodecHeader readHeader(PositionalBufferedStream stream) throws IOException {
        return this.wrappedCodec.readHeader(stream);
    }

    @Override
    public Class<VCFVariant> getFeatureType() {
        return VCFVariant.class;
    }

    @Override
    public boolean canDecode(String path) {
        return path.endsWith(".bcf");
    }
}
