package org.igv.hic;

import org.igv.ucsc.twobit.UnsignedByteBuffer;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * Static block index stored in memory.
 * Uses BinaryParser to read entries: int blockNumber, long filePosition, int size.
 */
public class StaticBlockIndex {

    private final Map<Integer, BlockIndexEntry> blockIndex = new HashMap<>();

    public StaticBlockIndex(int nBlocks, UnsignedByteBuffer dis) throws IOException {
        while (nBlocks-- > 0) {
            int blockNumber = dis.getInt();
            long filePosition = dis.getLong();
            int size = dis.getInt();
            blockIndex.put(blockNumber, new BlockIndexEntry(filePosition, size));
        }
    }

    public BlockIndexEntry getBlockIndexEntry(int blockNumber) {
        return blockIndex.get(blockNumber);
    }

    public static class BlockIndexEntry {
        public final long filePosition;
        public final int size;
        public BlockIndexEntry(long filePosition, int size) {
            this.filePosition = filePosition;
            this.size = size;
        }
    }
}