package org.broad.igv.feature.genome;

import java.io.IOException;
import java.util.List;

public class TwoBitSequence implements Sequence{

    private TwoBitReader reader;

    public TwoBitSequence(String path) throws IOException {
        this.reader = new TwoBitReader(path);
    }

    @Override
    public byte[] getSequence(String chr, int start, int end) {
        return this.reader.readSequence(chr, start, end);
    }

    @Override
    public byte getBase(String chr, int position) {
        throw new RuntimeException("getBase is not implementd for TwoBitSequence");
    }

    @Override
    public List<String> getChromosomeNames() {
         return this.reader.getSequenceNames();
    }

    @Override
    public int getChromosomeLength(String chrname) {
        throw new RuntimeException("getBase is not implementd for TwoBitSequence");
    }

    @Override
    public boolean isRemote() {
        return false;
    }
}
