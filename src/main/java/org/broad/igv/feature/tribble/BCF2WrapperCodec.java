package org.broad.igv.feature.tribble;

import htsjdk.samtools.util.LocationAware;
import org.broad.igv.logging.*;
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

    private static Logger log = LogManager.getLogger(BCF2WrapperCodec.class);

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
        String chr = genome == null ? vc.getChr() : genome.getCanonicalChrName(vc.getChr());
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
