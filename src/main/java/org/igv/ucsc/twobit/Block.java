package org.igv.ucsc.twobit;

class Block {
    final long start;
    final int size;
    final long end;

    Block(long start, int size) {
        this.start = start;
        this.size = size;
        this.end = start + size;
    }
}
