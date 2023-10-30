package org.broad.igv.ucsc;


import htsjdk.samtools.seekablestream.SeekableStream;
import org.broad.igv.util.UnsignedByteBuffer;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by jrobinso on 6/13/17.
 */
public class RPTree implements BPIndex{

    // the number 0x78CA8C91 in the architecture of the machine that created the file
    static int SIGNATURE = 0x78CA8C91;

    SeekableStream is;
    private long fileOffset;
    ByteOrder byteOrder = ByteOrder.LITTLE_ENDIAN;  // Until proven otherwise

    Map<Long, Node> nodeCache;

    // Header
    int blockSize;
    int keySize;
    int valSize;
    long itemCount;
    long reserved;
    long nodeOffset;

    public RPTree(String path, long fileOffset) throws IOException {
        this.is = IGVSeekableStreamFactory.getInstance().getStreamFor(path);
        this.fileOffset = fileOffset;
        this.nodeCache = new HashMap<>();
        init();
    }

    UnsignedByteBuffer loadBinaryBuffer(long start, int size) throws IOException {
        ByteBuffer bb = ByteBuffer.allocate(size);
        bb.order(this.byteOrder);
        byte[] bytes = bb.array();
        this.is.seek(start);
        this.is.readFully(bytes);
        return new UnsignedByteBuffer(bb);
    }

    private void init() throws IOException {

        long filePosition = this.fileOffset;
        UnsignedByteBuffer buffer = loadBinaryBuffer(filePosition, 64);

        int magicNumber = buffer.getInt();
        if (SIGNATURE != magicNumber) {
            this.byteOrder = ByteOrder.BIG_ENDIAN;
            buffer = loadBinaryBuffer(0, 64);  // Reload for new byte order
            magicNumber = buffer.getInt();
            if (SIGNATURE != magicNumber) {
                throw new RuntimeException("Unexpected magic number");
            }
        }
        this.blockSize = buffer.getInt();
        this.keySize = buffer.getInt();
        this.valSize = buffer.getInt();
        this.itemCount = buffer.getLong();
        this.reserved = buffer.getLong();
        this.nodeOffset = this.fileOffset + 32;
    }

    public long[] search(String term) throws IOException {

        int keySize = this.keySize;
        int valSize = this.valSize;

        if (!(valSize == 16 || valSize == 8)) {
            throw new RuntimeException("Unexpected valSiz: " + this.valSize);
        }

        // Kick things off
        return walkTreeNode(this.nodeOffset, term);
    }

    Node readTreeNode(long offset) throws IOException {

        if (this.nodeCache.containsKey(offset)) {
            return this.nodeCache.get(offset);
        } else {
            UnsignedByteBuffer buffer = loadBinaryBuffer(offset, 4);

            byte type = buffer.get();
            byte reserved = buffer.get();
            int count = buffer.getUShort();
            List<Item> items = new ArrayList<>();

            if (type == 1) {
                // Leaf node
                int size = count * (keySize + valSize);

                buffer = loadBinaryBuffer(offset + 4, size);
                for (int i = 0; i < count; i++) {
                    String key = buffer.getFixedLengthString(keySize);
                    long itemOffset = buffer.getLong();
                    long[] value;
                    if (valSize == 16) {
                        long length = buffer.getLong();
                        value = new long[]{itemOffset, length};
                    } else {
                        value = new long[]{itemOffset};
                    }
                    items.add(new Item(key, value));
                }
            } else {
                // Non leaf node
                int size = count * (keySize + 8);
                buffer = loadBinaryBuffer(offset + 4, size);

                for (int i = 0; i < count; i++) {
                    String key = buffer.getFixedLengthString(keySize);
                    long itemOffset = buffer.getLong();
                    long[] value = new long[]{itemOffset};
                    items.add(new Item(key, value));
                }
            }

            Node node = new Node(type, count, items);
            this.nodeCache.put(offset, node);
            return node;
        }
    }

    long[] walkTreeNode(long offset, String term) throws IOException {

        Node node = readTreeNode(offset);

        if (node.type == 1) {
            // Leaf node
            for (Item item : node.items) {
                if (term.equals(item.key)) {
                    return item.value;
                }
            }
            // If we get here, not found
            return null;
        } else {

            // Non leaf node

            // Read and discard the first key.
            long childOffset = node.items.get(0).value[0];

            for (int i = 1; i < node.items.size(); i++) {
                String key = node.items.get(i).key;
                if (term.compareTo(key) < 0) {
                    break;
                }
                childOffset = node.items.get(i).value[0];
            }
            return walkTreeNode(childOffset, term);
        }
    }


    class Node {
        int type;
        int count;
        List<Item> items;

        public Node(int type, int count, List<Item> items) {
            this.type = type;
            this.count = count;
            this.items = items;
        }
    }

    class Item {
        String key;
        long[] value;

        public Item(String key, long[] value) {
            this.key = key;
            this.value = value;
        }
    }

}






