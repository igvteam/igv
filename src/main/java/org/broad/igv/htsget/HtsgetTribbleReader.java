package org.broad.igv.htsget;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.FeatureCodec;

import java.io.IOException;
import java.util.List;

public class HtsgetTribbleReader extends AbstractFeatureReader {

    public HtsgetTribbleReader(String path, FeatureCodec codec) {
        super(path, codec);
    }

    @Override
    public CloseableTribbleIterator query(String chr, int start, int end) throws IOException {
        return null;
    }

    @Override
    public CloseableTribbleIterator iterator() throws IOException {
        return null;
    }

    @Override
    public void close() throws IOException {

    }

    @Override
    public List<String> getSequenceNames() {
        return null;
    }
}
