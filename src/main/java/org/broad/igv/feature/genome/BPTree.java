package org.broad.igv.feature.genome;

/**
 * A UCSC (Jim Kent) B+ Tree.
 * From Kent  https://github.com/ucscGenomeBrowser/kent/blob/master/src/lib/bPlusTree.c
 * *
 * * The layout of the file on disk is:
 * *    header
 * *    root node
 * *    (other nodes)
 * * In general when the tree is first built the higher level nodes are stored before the
 * * lower level nodes.  It is possible if a b+ tree is dynamically updated for this to
 * * no longer be strictly true, but actually currently the b+ tree code here doesn't implement
 * * dynamic updates - it just creates a b+ tree from a sorted list.
 * *
 * * Each node can be one of two types - index or leaf.  The index nodes contain pointers
 * * to child nodes.  The leaf nodes contain the actual data.
 * *
 * * The layout of the file header is:
 * *       <magic number>  4 bytes - The value bptSig (0x78CA8C91)
 * *       <block size>    4 bytes - Number of children per block (not byte size of block)
 * *       <key size>      4 bytes - Number of significant bytes in key
 * *       <val size>      4 bytes - Number of bytes in value
 * *       <item count>    8 bytes - Number of items in index
 * *       <reserved2>     4 bytes - Always 0 for now
 * *       <reserved3>     4 bytes - Always 0 for now
 * * The magic number may be byte-swapped, in which case all numbers in the file
 * * need to be byte-swapped.
 * *
 * * The nodes start with a header:
 * *       <is leaf>       1 byte  - 1 for leaf nodes, 0 for index nodes.
 * *       <reserved>      1 byte  - Always 0 for now.
 * *       <count>         2 bytes - The number of children/items in node
 * * This is followed by count items.  For the index nodes the items are
 * *       <key>           key size bytes - always written most significant byte first
 * *       <offset>        8 bytes - Offset of child node in index file.
 * * For leaf nodes the items are
 * *       <key>           key size bytes - always written most significant byte first
 * *       <val>           val sized bytes - the value associated with the key.
 * * Note in general the leaf nodes may not be the same size as the index nodes, though in
 * * the important case where the values are file offsets they will be.
 */


import htsjdk.samtools.seekablestream.SeekableStream;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.*;

/**
 * Created by jrobinso on 6/13/17.
 */
public class BPTree {

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

    public BPTree(String path, long fileOffset) throws IOException {
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

    long[] search(String term) throws IOException {

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






