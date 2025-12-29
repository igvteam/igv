package org.igv.ucsc.twobit;

import java.io.IOException;

public class SequenceRecord {

    final int dnaSize;
    Block[] nBlocks;
    Block[] maskBlocks;

    final long packedPos;

    SequenceRecord(int dnaSize, Block[] nBlocks, Block[] maskBlocks, long packedPos) throws IOException {
        this.dnaSize = dnaSize;
        this.nBlocks = nBlocks;
        this.maskBlocks = maskBlocks;
        this.packedPos = packedPos;
    }

    public int getDnaSize() {
        return dnaSize;
    }
}
