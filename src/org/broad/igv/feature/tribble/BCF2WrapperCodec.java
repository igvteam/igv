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

import htsjdk.samtools.util.LocationAware;
import org.apache.log4j.Logger;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.variant.vcf.VCFVariant;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.variantcontext.VariantContext;

import java.io.IOException;
import java.io.InputStream;

/**
 * @author Jacob Silterra
 * @date 2013-Jun-14
 */
public class BCF2WrapperCodec implements FeatureCodec<VCFVariant, PositionalBufferedStream> {

    private static Logger log = Logger.getLogger(BCF2WrapperCodec.class);

    FeatureCodec<VariantContext, PositionalBufferedStream> wrappedCodec;
    Genome genome;

    public BCF2WrapperCodec(FeatureCodec<VariantContext, PositionalBufferedStream> wrappedCodec, Genome genome) {
        this.wrappedCodec = wrappedCodec;
        this.genome = genome;
    }

    @Override
    public VCFVariant decode(PositionalBufferedStream stream) throws IOException {
        VariantContext vc = wrappedCodec.decode(stream);
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

    @Override
    public LocationAware makeIndexableSourceFromStream(final InputStream bufferedInputStream) {
        final PositionalBufferedStream pbs;
        if (bufferedInputStream instanceof PositionalBufferedStream) {
            pbs = (PositionalBufferedStream) bufferedInputStream;
        } else {
            pbs = new PositionalBufferedStream(bufferedInputStream);
        }
        return new AsciiLineReaderIterator(new AsciiLineReader(pbs));
    }

    @Override
    public boolean isDone(PositionalBufferedStream positionalBufferedStream) {
        try {
            return positionalBufferedStream.isDone();
        } catch (IOException e) {
            log.error(e.getMessage(), e);
            return true;
        }
    }

    @Override
    public void close(PositionalBufferedStream positionalBufferedStream) {
        positionalBufferedStream.close();
    }

    @Override
    public PositionalBufferedStream makeSourceFromStream(final InputStream bufferedInputStream) {
        return new PositionalBufferedStream(bufferedInputStream);
    }
}
