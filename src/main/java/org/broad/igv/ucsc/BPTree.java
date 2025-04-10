package org.broad.igv.ucsc;

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


import org.broad.igv.ucsc.twobit.TwoBitSequence;
import org.broad.igv.ucsc.twobit.UnsignedByteBuffer;
import org.broad.igv.ucsc.twobit.UnsignedByteBufferImpl;

import java.io.IOException;
import java.nio.ByteOrder;
import java.util.*;

/**
 * Created by jrobinso on 6/13/17.
 */
public class BPTree implements BPIndex{

    // the number 0x78CA8C91 in the architecture of the machine that created the file
    static int SIGNATURE = 0x78CA8C91;

    String path;
    private long fileOffset;
    ByteOrder byteOrder = ByteOrder.LITTLE_ENDIAN;  // Until proven otherwise

    Map<Long, Node> nodeCache;

    // Header
    int blockSize;
    int keySize;
    int valSize;
    public long itemCount;
    long reserved;
    long nodeOffset;

    public static BPTree loadBPTree(String path, long fileOffset) throws IOException {
        BPTree tree = new BPTree(path, fileOffset);
        tree.init();
        return tree;
    }

    private BPTree(String path, long fileOffset) throws IOException {
        this.path = path;
        this.fileOffset = fileOffset;
        this.nodeCache = new HashMap<>();
    }

    UnsignedByteBuffer loadBinaryBuffer(long start, int size) throws IOException {
        return UnsignedByteBufferImpl.loadBinaryBuffer(this.path, this.byteOrder, start, size);
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
        if (!(valSize == 16 || valSize == 8)) {
            throw new RuntimeException("Unexpected valSiz: " + this.valSize);
        }
        this.itemCount = buffer.getLong();
        this.reserved = buffer.getLong();
        this.nodeOffset = this.fileOffset + 32;
    }

    private byte[] search(String term) throws IOException {

        // Kick things off
        return walkTreeNode(this.nodeOffset, term);
    }

    public long searchForOffset(String term) throws IOException {
        byte[] bytes = search(term);
        if (bytes != null) {
            return bytesToLong(bytes, 0);
        } else {
            return -1;
        }
    }

    public long[] searchLongLong(String term) throws IOException {
        byte[] bytes = search(term);
        if (bytes != null) {
            long offset = bytesToLong(bytes, 0);
            long size = bytesToLong(bytes, 8);
            return new long[]{offset, size};
        } else {
            return null;
        }
    }

    public int []  searchIntInt(String term) throws IOException {
        byte[] bytes = search(term);
        if (bytes != null) {
            int offset = bytesToInt(bytes, 0);
            int size = bytesToInt(bytes, 4);
            return new int[]{offset, size};
        } else {
            return null;
        }
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
                    byte[] value = buffer.getBytes(valSize);
                    items.add(new Item(key, value));
                }
            } else {
                // Non leaf node
                int size = count * (keySize + 8);
                buffer = loadBinaryBuffer(offset + 4, size);

                for (int i = 0; i < count; i++) {
                    String key = buffer.getFixedLengthString(keySize);
                    long childOffset = buffer.getLong();
                    items.add(new Item(key, childOffset));
                }
            }

            Node node = new Node(type, count, items);
            this.nodeCache.put(offset, node);
            return node;
        }
    }

    byte[] walkTreeNode(long offset, String term) throws IOException {

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
            long childOffset = node.items.get(0).offset;
            for (int i = 1; i < node.items.size(); i++) {
                String key = node.items.get(i).key;
                if (term.compareTo(key) < 0) {
                    break;
                }
                childOffset = node.items.get(i).offset;
            }
            return walkTreeNode(childOffset, term);
        }
    }


    private long bytesToLong(byte[] bytes, int start) {
        long value = 0;
        int length = Math.min(bytes.length, 8);
        if (byteOrder == ByteOrder.BIG_ENDIAN) {
            for (int i = start; i < start + length; i++) {
                value <<= 8;
                value |= (bytes[i] & 0xFF);
            }
        } else { // Little-endian
            for (int i = start + length - 1; i >= start; i--) {
                value <<= 8;
                value |= (bytes[i] & 0xFF);
            }
        }
        return value;
    }

    private int bytesToInt(byte[] bytes, int start) {
        int value = 0;
        if (byteOrder == ByteOrder.BIG_ENDIAN) {
            for (int i = start; i < start + 4; i++) {
                value <<= 8;
                value |= (bytes[i] & 0xFF);
            }
        } else { // Little-endian
            for (int i = start + 4 - 1; i >= start; i--) {
                value <<= 8;
                value |= (bytes[i] & 0xFF);
            }
        }
        return value;
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

    /**
     * Represents a single item in a B+ tree node. An item should have a value (for leaf nodes),  or an offset, but not both.
     */
    class Item {
        String key;
        byte[] value;
        long offset;

        public Item(String key, long offset) {
            this.key = key;
            this.offset = offset;
        }

        public Item(String key, byte[] value) {
            this.key = key;
            this.value = value;
        }
    }
}







