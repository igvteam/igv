package org.broad.igv.ucsc.bb;


import org.broad.igv.ucsc.twobit.UnsignedByteBuffer;
import org.broad.igv.ucsc.twobit.UnsignedByteBufferImpl;

import java.io.IOException;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class RPTree {
    static int RPTREE_HEADER_SIZE = 48;
    static int RPTREE_NODE_LEAF_ITEM_SIZE = 32; // leaf item size
    static int RPTREE_NODE_CHILD_ITEM_SIZE = 24; // child item size

    static int magic = 610839776;

    ByteOrder byteOrder = ByteOrder.LITTLE_ENDIAN;

    Map<Long, Node> nodeCache = new HashMap<>();
    long startOffset;
    String path;
    private Header header;
    private long rootNodeOffset;

    static RPTree loadTree(String path, long startOffset) throws IOException {
        RPTree tree = new RPTree(path, startOffset);
        tree.init();
        return tree;
    }
    private RPTree(String path, long startOffset) {
        this.path = path;
        this.startOffset = startOffset;
    }

    void init() throws IOException {
        UnsignedByteBuffer binaryParser = UnsignedByteBufferImpl.loadBinaryBuffer(this.path, this.byteOrder, this.startOffset, RPTREE_HEADER_SIZE);
        int magic = binaryParser.getInt();
        if (magic != RPTree.magic) {
            this.byteOrder = ByteOrder.BIG_ENDIAN;
            binaryParser = UnsignedByteBufferImpl.loadBinaryBuffer(this.path, this.byteOrder, this.startOffset, RPTREE_HEADER_SIZE);
            magic = binaryParser.getInt();
            if (magic != RPTree.magic) {
                throw new RuntimeException("Bad magic number " + magic);
            }
        }
        this.header = new Header(binaryParser);
        this.rootNodeOffset = this.startOffset + RPTREE_HEADER_SIZE;
    }

/**
     * Find all leaf items in the tree that overlap the specified region.
     *
     * @param chrIdx1 The first chromosome index
     * @param startBase The start base of the region
     * @param chrIdx2 The second chromosome index
     * @param endBase The end base of the region
     * @return A list of leaf items that overlap the specified region
     */
    List<Item> findLeafItemsOverlapping(int chrIdx1, int startBase, int chrIdx2, int endBase) throws IOException {
        List<Item> leafItems = new ArrayList<>();
        this.walkTree(this.rootNodeOffset, leafItems, chrIdx1, startBase, chrIdx2, endBase);
        return leafItems;
    }

    /**
     * Find all leaf items in the tree.
     *
     * @return A list of all leaf items
     */
    List<Item> findAllLeafItems() throws IOException {
        List<Item> leafItems = new ArrayList<>();
        this.walkTree(this.rootNodeOffset, leafItems, -1, -1, -1, -1);
        return leafItems;
    }

    /**
     * Walk the tree recursively, adding items to the list if they overlap the specified region.
     *
     * @param offset      The offset of the current node
     * @param leafItems   The list of items to add to
     * @param chrIdx1     The first chromosome index
     * @param startBase   The start base of the region
     * @param chrIdx2     The second chromosome index
     * @param endBase     The end base of the region
     */
    void walkTree(long offset, List<Item> leafItems, int chrIdx1, int startBase, int chrIdx2, int endBase) throws IOException {
            Node node = readNode(offset);
            for (Item item : node.items) {
                if (chrIdx1 < 0 || item.overlaps(chrIdx1, startBase, chrIdx2, endBase)) {
                    if (node.type == 1) {   // Leaf node
                        leafItems.add(item);
                    } else { // Non leaf node
                        this.walkTree(item.dataOffset, leafItems, chrIdx1, startBase, chrIdx2, endBase);
                    }
                }
            }
    }

    static class Header {
        long blockSize,
                itemCount,
                startChromIx,
                startBase,
                endChromIx,
                endBase,
                endFileOffset,
                itemsPerSlot,
                reserved;

        Header(UnsignedByteBuffer binaryParser) {
            blockSize = binaryParser.getUInt();
            itemCount = binaryParser.getLong();
            startChromIx = binaryParser.getUInt();
            startBase = binaryParser.getUInt();
            endChromIx = binaryParser.getUInt();
            endBase = binaryParser.getUInt();
            endFileOffset = binaryParser.getLong();
            itemsPerSlot = binaryParser.getUInt();
            reserved = binaryParser.getUInt();
        }
    }

    Node readNode(long offset) throws IOException {

        long nodeKey = offset;
        if (this.nodeCache.containsKey(nodeKey)) {
            return this.nodeCache.get(nodeKey);
        }

        UnsignedByteBuffer binaryParser = UnsignedByteBufferImpl.loadBinaryBuffer(this.path, this.byteOrder, offset, 4);
        byte type = binaryParser.get();
        boolean isLeaf = (type == 1);
        byte reserved = binaryParser.get();
        int count = binaryParser.getUShort();
        int bytesRequired = count * (isLeaf ? RPTREE_NODE_LEAF_ITEM_SIZE : RPTREE_NODE_CHILD_ITEM_SIZE);
        binaryParser = UnsignedByteBufferImpl.loadBinaryBuffer(this.path, this.byteOrder, offset + 4, bytesRequired);

        Item[] items = new Item[count];
        for (int i = 0; i < count; i++) {
            items[i] = new Item(binaryParser, type);
        }
        Node node = new Node(type, items);
        this.nodeCache.put(nodeKey, node);
        return node;
    }

    static class Node {
        int type;
        Item[] items;

        public Node(int type, Item[] items) {
            this.type = type;
            this.items = items;
        }
    }

    static class Item {
        int startChrom;
        int startBase;
        int endChrom;
        int endBase;
        long dataOffset;
        long dataSize;

        Item(UnsignedByteBuffer binaryParser, int type) {

            startChrom = binaryParser.getInt();
            startBase = binaryParser.getInt();
            endChrom = binaryParser.getInt();
            endBase = binaryParser.getInt();
            dataOffset = binaryParser.getLong();
            if (type == 1) {
                dataSize = binaryParser.getLong();
            }

        }

        /**
         * Return true if {chrIdx1:startBase-chrIdx2:endBase} overlaps item's interval
         *
         * @returns {boolean}
         */
        boolean overlaps(int chrIdx1, int startBase, int chrIdx2, int endBase) {

            return ((chrIdx2 > this.startChrom) || (chrIdx2 == this.startChrom && endBase >= this.startBase)) &&
                    ((chrIdx1 < this.endChrom) || (chrIdx1 == this.endChrom && startBase <= this.endBase));


        }
    }


}

