package org.igv.hic;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.util.HashMap;
import java.util.Map;

/**
 * Dynamic block index that lazily reads index entries from a file.
 * Constructor expects a FileChannel positioned at block index start.
 */
public class DynamicBlockIndex {

    private final FileChannel fileChannel;
    private final long minPosition;
    private final long maxPosition;
    private final int nBlocks;
    private final int maxBlock;

    private Map<Integer, StaticBlockIndex.BlockIndexEntry> blockIndexMap;
    private Range blockNumberRange;
    private MapFileBounds mapFileBounds;

    private static final long ENTRY_SIZE = 16L;

    public DynamicBlockIndex(FileChannel fileChannel, long position, int nBlocks, int maxBlock) {
        this.fileChannel = fileChannel;
        this.minPosition = position;
        this.maxPosition = position + (long) nBlocks * ENTRY_SIZE;
        this.nBlocks = nBlocks;
        this.maxBlock = maxBlock;
    }

    public StaticBlockIndex.BlockIndexEntry getBlockIndexEntry(int blockNumber) throws IOException {
        if (blockNumber > maxBlock) {
            return null;
        } else if (blockNumberRange != null && blockNumber >= blockNumberRange.first && blockNumber <= blockNumberRange.last) {
            return blockIndexMap.get(blockNumber);
        } else {
            long minPos = minPosition;
            long maxPos = maxPosition;
            if (blockNumberRange != null && mapFileBounds != null) {
                if (blockNumber < blockNumberRange.first) {
                    maxPos = mapFileBounds.min;
                } else if (blockNumber > blockNumberRange.last) {
                    minPos = mapFileBounds.max;
                }
            }
            if (maxPos - minPos < ENTRY_SIZE) {
                return null;
            } else {
                return searchForBlockIndexEntry(blockNumber, minPos, maxPos);
            }
        }
    }

    private StaticBlockIndex.BlockIndexEntry searchForBlockIndexEntry(int blockNumber, long boundsMin, long boundsMax) throws IOException {
        final long chunkSize = ENTRY_SIZE * 100000L;
        if (boundsMax - boundsMin < chunkSize) {
            int len = (int) (boundsMax - boundsMin);
            ByteBuffer buf = ByteBuffer.allocate(len);
            buf.order(ByteOrder.LITTLE_ENDIAN);
            fileChannel.read(buf, boundsMin);
            buf.flip();

            Map<Integer, StaticBlockIndex.BlockIndexEntry> localIndex = new HashMap<>();
            int ptr = 0;
            Integer firstBlockNumber = null;
            Integer lastBlockNumber = null;
            while (ptr < len) {
                int bn = buf.getInt();
                long filePosition = buf.getLong();
                int blockSizeInBytes = buf.getInt();
                localIndex.put(bn, new StaticBlockIndex.BlockIndexEntry(filePosition, blockSizeInBytes));
                if (firstBlockNumber == null) firstBlockNumber = bn;
                lastBlockNumber = bn;
                ptr += ENTRY_SIZE;
            }
            this.mapFileBounds = new MapFileBounds(boundsMin, boundsMax);
            this.blockNumberRange = new Range(firstBlockNumber, lastBlockNumber);
            this.blockIndexMap = localIndex;
            return localIndex.get(blockNumber);
        }
        long nEntries = (boundsMax - boundsMin) / ENTRY_SIZE;
        long pos1 = boundsMin + (nEntries / 2) * ENTRY_SIZE;
        ByteBuffer buf = ByteBuffer.allocate((int) ENTRY_SIZE);
        buf.order(ByteOrder.LITTLE_ENDIAN);
        fileChannel.read(buf, pos1);
        buf.flip();
        int bn = buf.getInt();
        if (bn == blockNumber) {
            long filePosition = buf.getLong();
            int blockSizeInBytes = buf.getInt();
            return new StaticBlockIndex.BlockIndexEntry(filePosition, blockSizeInBytes);
        } else if (blockNumber > bn) {
            return searchForBlockIndexEntry(blockNumber, pos1 + ENTRY_SIZE, boundsMax);
        } else {
            return searchForBlockIndexEntry(blockNumber, boundsMin, pos1);
        }
    }

    private static class Range {
        final int first;
        final int last;

        Range(int first, int last) {
            this.first = first;
            this.last = last;
        }
    }

    private static class MapFileBounds {
        final long min;
        final long max;

        MapFileBounds(long min, long max) {
            this.min = min;
            this.max = max;
        }
    }
}