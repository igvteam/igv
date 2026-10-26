package org.igv.feature.tribble;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.variantcontext.VariantContext;
import org.igv.variant.vcf.VCFVariant;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;

/**
 * Reads a BCF file, raw or BGZF, as VCFVariants.  htsjdk queries a BGZF BCF through its CSI index only when given its
 * own BCF2Codec, which yields VariantContexts, so the conversion {@link BCF2WrapperCodec} makes per record for other
 * readers is made here, around the reader htsjdk chooses.
 */
public class VCFVariantReader extends AbstractFeatureReader<VCFVariant, PositionalBufferedStream> {

    private final AbstractFeatureReader<VariantContext, PositionalBufferedStream> reader;
    private final BCF2WrapperCodec codec;

    /**
     * @param indexPath a CSI index for BGZF BCF or a Tribble index for raw BCF, or null to look beside the file
     */
    public static VCFVariantReader open(String path, String indexPath, BCF2WrapperCodec codec, boolean requireIndex) {
        return new VCFVariantReader(path,
                AbstractFeatureReader.getFeatureReader(path, indexPath, codec.wrappedCodec, requireIndex), codec);
    }

    private VCFVariantReader(String path, AbstractFeatureReader<VariantContext, PositionalBufferedStream> reader,
                             BCF2WrapperCodec codec) {
        super(path, codec);
        this.reader = reader;
        this.codec = codec;
    }

    @Override
    public CloseableTribbleIterator<VCFVariant> query(String chr, int start, int end) throws IOException {
        return new VariantIterator(reader.query(chr, start, end));
    }

    @Override
    public CloseableTribbleIterator<VCFVariant> iterator() throws IOException {
        return new VariantIterator(reader.iterator());
    }

    @Override
    public List<String> getSequenceNames() {
        return reader.getSequenceNames();
    }

    @Override
    public Object getHeader() {
        return reader.getHeader();
    }

    @Override
    public boolean hasIndex() {
        return reader.hasIndex();
    }

    @Override
    public boolean isQueryable() {
        return reader.isQueryable();
    }

    @Override
    public void close() throws IOException {
        reader.close();
    }

    private class VariantIterator implements CloseableTribbleIterator<VCFVariant> {

        private final CloseableTribbleIterator<VariantContext> iterator;

        VariantIterator(CloseableTribbleIterator<VariantContext> iterator) {
            this.iterator = iterator;
        }

        @Override
        public boolean hasNext() {
            return iterator.hasNext();
        }

        @Override
        public VCFVariant next() {
            return codec.toVariant(iterator.next());
        }

        @Override
        public void close() {
            iterator.close();
        }

        @Override
        public Iterator<VCFVariant> iterator() {
            return this;
        }
    }
}
