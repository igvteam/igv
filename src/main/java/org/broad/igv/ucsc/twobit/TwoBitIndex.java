package org.broad.igv.ucsc.twobit;

import htsjdk.samtools.seekablestream.SeekableStream;
import org.broad.igv.ucsc.BPIndex;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

public class TwoBitIndex implements BPIndex {

    String path;
    ByteOrder byteOrder;

    Map<String, Long> sequenceDataOffsets;

    public TwoBitIndex(String path, ByteOrder byteOrder, int seqCount) throws IOException {
        this.path = path;
        this.byteOrder = byteOrder;
        this.sequenceDataOffsets = new HashMap<>();
        readIndex(seqCount);
    }

    private void readIndex(int seqCount) throws IOException {

        long filePosition = 16;

        int estNameLength = 20;
        int estSize = seqCount * estNameLength + 100;
        ByteBuffer buffer = loadBinaryBuffer(filePosition, estSize);

        // Loop through sequences loading name and file offset.  We don't know the precise size in bytes in advance
        // so we need to check for bytes available and reload as needed.

        sequenceDataOffsets = new LinkedHashMap<>();
        for (int i = 0; i < seqCount; i++) {

            if (buffer.remaining() < 1) {
                filePosition += buffer.position();
                estSize = (seqCount - i) * estNameLength + 100;
                buffer = loadBinaryBuffer(filePosition, estSize);
            }

            final byte nameSize = buffer.get();

            if (buffer.remaining() < nameSize * 5) {
                filePosition += buffer.position();
                estSize = (seqCount - i) * estNameLength + 100;
                buffer = loadBinaryBuffer(filePosition, estSize);
            }

            byte[] seqNameBytes = new byte[nameSize];
            buffer.get(seqNameBytes);
            String seqName = new String(seqNameBytes);

            long offset = buffer.getInt();
            sequenceDataOffsets.put(seqName, offset);
        }
    }

    public long[] search(String term) {
        if(sequenceDataOffsets.containsKey(term)) {
            return new long[]{(long) sequenceDataOffsets.get(term)};
        } else {
            return null;
        }
    }

    ByteBuffer loadBinaryBuffer(long start, int size) throws IOException {
        try (SeekableStream is = IGVSeekableStreamFactory.getInstance().getStreamFor(path)) {
            ByteBuffer bb = ByteBuffer.allocate(size);
            bb.order(this.byteOrder);
            byte[] bytes = bb.array();
            is.seek(start);
            is.readFully(bytes);
            return bb;
        }
    }
}
