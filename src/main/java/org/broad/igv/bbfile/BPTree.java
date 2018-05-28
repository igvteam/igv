/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.bbfile;

import htsjdk.samtools.seekablestream.SeekableStream;
import org.apache.log4j.Logger;
import org.broad.igv.util.LittleEndianInputStream;

import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: martind
 * Date: Dec 17, 2009
 * Time: 12:28:30 PM
 * To change this template use File | Settings | File Templates.
 */
/*
*   B+ Tree class will construct a B+ tree from a binary Bed/Wig BBFile.
*   (or by insertion of tree nodes - TBD see insert method)
*
*   1) BPTree will first read in the B+ tree header with BPTreeHeader class.
*
*   2) Starting with the root node, the readBPTreeNode method will read in the
*   node format, determine if the node contains child nodes (isLeaf = false)
*   or leaf items (isLeaf = true).
*
*   3) If node is a leaf node, all leaf items are read in to the node's leaf array.
*
*   4) If node is a child node, readBPTreeNode will be called recursively,
*   until the leaf node is encountered, where step 3 is performed.
*
*   5) The child nodes will be populated with their child node items in reverse order
*   of recursion from step 4, until the tree is completely populated
*   back up to the root node.
*
*   6) The getChromosomeKey is provided to construct a valid key for B+
*   chromosome tree searches, and getChromosomeID returns a chromosome ID for
*   searches in the R+ index tree.
*
**/
public class BPTree {

    private static Logger log = Logger.getLogger(BPTree.class);

    public static final int BPTREE_NODE_FORMAT_SIZE = 4;   // node format size
    public static final int BPTREE_NODE_ITEM_SIZE = 8;     // Plus keySize to be added

    // B+ tree access variables   - for reading in B+ tree nodes from a file
    private SeekableStream fis;      // file handle - BBFile input stream
    private long treeOffset;         // mChromosome B+ tree file offset
    private BPTreeHeader treeHeader; // B+ tree header (Table E for BBFile)

    // B+ tree organizational variables  - derived from Table E
    private int blockSize;     // number of children per block
    private int keySize;       // character size of primary key
    private int valueSize;     // number of bytes in value being indexed
    private long itemCount;    //  number of contig/mChromosome items in tree

    // B+ tree nodal variables
    private BPTreeNode rootNode;   // B+ tree root node
    private long nodeCount;        // number of nodes defined in the B+ tree
    private long leafCount;        // number of leaves in the B+ tree

    private Map<Integer, String> idChromMap;
    private Map<String, Integer> chromIdMap;


    /*
    *    Constructor for reading in a B+ tree from a BBFile/input stream.
    *
    *    Parameters:
    *        fis - file input stream handle
    *        fileOffset - file offset to the B+ tree header
    *        isLowToHigh - indicates byte order is low to high, else is high to low
    * */
    public BPTree(SeekableStream fis, long fileOffset, boolean isLowToHigh) {

        // Save the seekable file handle and B+ Tree file offset
        // Note: the offset is the B+ Tree Header Table E file location
        this.fis = fis;
        treeOffset = fileOffset;

        idChromMap = new HashMap<>();
        chromIdMap = new HashMap<>();

        // read in B+ tree header - verify the B+ tree info exits
        treeHeader = new BPTreeHeader(this.fis, treeOffset, isLowToHigh);

        // log error if header not found and throw exception
        if (!treeHeader.isHeaderOK()) {
            int badMagic = treeHeader.getMagic();
            log.error("Error reading B+ tree header: bad magic = " + badMagic);
            throw new RuntimeException("Error reading B+ tree header: bad magic = "
                    + badMagic);
        }

        // assign B+ tree specifications from the header
        blockSize = treeHeader.getBlockSize();
        keySize = treeHeader.getKeySize();
        valueSize = treeHeader.getValSize();
        itemCount = treeHeader.getItemCount();

        // populate the tree - read in the nodes
        long nodeOffset = treeOffset + treeHeader.BPTREE_HEADER_SIZE;
        BPTreeNode parentNode = null;  // parent node of the root is itself, or null

        // get the root node - which recursively populates the remaining nodes
        rootNode = readBPTreeNode(this.fis, nodeOffset, parentNode, isLowToHigh);

    }

    /*
    *   Method returns the file input stream handle
    * */
    public SeekableStream getFis() {
        return fis;
    }

    /*
    *   Method returns the chromosome name key size, which is
    *   the number of valid characters for chromosome name.
    * */
    public int getKeySize() {
        return keySize;
    }

    /*
    * Method returns the number of chromosome/contig names.
    */
    public long getItemCount() {
        return itemCount;
    }


    Map<String, String> chromosomeKeyCache = new HashMap();

    /*
    *   Returns a search key for the mChromosome region  which  can
    *   be used to search for a corresponding section in the B+ tree.
    *
    *   According the the spec the key is the "first keySize characters of chromosome name"
    * */
    public String getChromosomeKey(String chromosome) {

        String key = chromosomeKeyCache.get(chromosome);
        if (key == null) {
            key = chromosome.length() <= keySize ? chromosome : chromosome.substring(0, keySize);
            chromosomeKeyCache.put(chromosome, key);
        }
        return key;
    }

    /*
    *   Returns a chromosome ID  which  can be used to search for a
    *   corresponding data section in the R+ tree for data.
    *
       Parameters:
    *       chromKey - chromosome name of valid key size.
    *
    *
    *   Note: A chromosomeID of -1 means chromosome name not included in B+ tree.
    *
    * */
    public int getChromosomeID(String chromKey) {
        return chromIdMap.get(chromKey);
    }

    /*
    *   Returns a chromosome name which is the B+ key for returning the
    *   chromosome ID for lookup in the R+ tree for data.
    *
    *   Parameters:
    *       chromID - chromosome ID expected in B+ tree
    *
    *   Returns:
    *       Chromosome name key; a null string means chromosome ID not found.
    *
    * */
    public String getChromosomeName(int chromID) {

        return idChromMap.get(chromID);
    }

    /*
    *   Method returns all chromosome key names in B+ tree.
    *
    *   Returns:
    *   Collection of all (chromosome ID, chromosome name)entries
    * */
    public ArrayList<String> getChromosomeNames() {

        return new ArrayList<>(chromIdMap.keySet());
    }

    /*
   *   Method returns all chromosome name, chromosome ID pairs for a given ID range.
   *
   *   Parameters:
   *       startChromID - starting ID for chromosome range expected in B+ tree
   *       endChromID - ending ID for chromosome range expected in B+ tree
   *
   *   Returns:
   *       Collection of (chromosome ID, chromosome name key) hash items;
   *       where an empty collection means ID range was not found.
   *
   * */
    public Map<Integer, String> getChromosomeIDMap(int startChromID, int endChromID) {

        Map<Integer, String> m = new HashMap<>();
        for (int i = startChromID; i <= endChromID; i++) {
            String c = idChromMap.get(i);
            if (c != null) {
                m.put(i, c);
            }
        }
        return m;
    }

    // prints out the B+ Tree  nodes and leaves
    public void print() {

        // check if read in
        if (!treeHeader.isHeaderOK()) {
            int badMagic = treeHeader.getMagic();
            log.error("Error reading B+ tree header: bad magic = " + badMagic);
            return;
        }

        // print B+ tree header
        treeHeader.print();

        // print  B+ tree node and leaf items - recursively
        if (rootNode != null)
            rootNode.printItems();
    }

    /*
    *   Method reads in the B+ tree nodes from the file, recursively.
    *
    *   Parameters:
    *       fis - file input stream handle
    *       fileOffset - file offset for B+ tree header
    *       keySize - chromosome name key size in characters
    *       parent - parent node
    *       isLowToHigh - if true, indicates byte order is low to high; else is high to low
    *
    *   Returns:
     *      Boolean which indicates if the B+ tree header was read correctly, with
    *       true for success, false for failure to find the header information.
    * */
    private BPTreeNode readBPTreeNode(SeekableStream fis, long fileOffset,
                                      BPTreeNode parent, boolean isLowToHigh) {

        LittleEndianInputStream lbdis = null;     // low to high byte reader
        DataInputStream bdis = null;        // high to low byte reader

        // set up for node format
        byte[] buffer = new byte[BPTREE_NODE_FORMAT_SIZE];
        BPTreeNode thisNode = null;
        BPTreeNode childNode = null;

        byte type;
        byte bval;
        int itemCount;
        int itemSize;
        boolean isLeaf;

        try {

            // Read node format into a buffer
            fis.seek(fileOffset);
            fis.readFully(buffer);

            if (isLowToHigh)
                lbdis = new LittleEndianInputStream(new ByteArrayInputStream(buffer));
            else
                bdis = new DataInputStream(new ByteArrayInputStream(buffer));

            // find node type
            if (isLowToHigh)
                type = lbdis.readByte();
            else
                type = bdis.readByte();

            // create the B+ tree node
            if (type == 1) {
                isLeaf = true;
                thisNode = new BPTreeLeafNode(++nodeCount);
            } else {
                isLeaf = false;
                thisNode = new BPTreeChildNode(++nodeCount);
            }

            if (isLowToHigh) {
                bval = lbdis.readByte();      // reserved - not currently used
                itemCount = lbdis.readUShort();
            } else {
                bval = bdis.readByte();      // reserved - not currently used
                itemCount = bdis.readUnsignedShort();
            }

            // Note: B+ tree node item size is the same for leaf and child items
            itemSize = BPTREE_NODE_ITEM_SIZE + this.keySize;
            int totalSize = itemSize * itemCount;
            byte[] itemBuffer = new byte[totalSize];
            fis.readFully(itemBuffer);

            if (isLowToHigh)
                lbdis = new LittleEndianInputStream(new ByteArrayInputStream(itemBuffer));
            else
                bdis = new DataInputStream(new ByteArrayInputStream(itemBuffer));

            // get the node items - leaves or child nodes
            for (int item = 0; item < itemCount; ++item) {

                // always extract the key from the node format
                char[] keychars = new char[keySize];  // + 1 for 0 byte
                int index;
                for (index = 0; index < keySize; ++index) {

                    if (isLowToHigh)
                        bval = lbdis.readByte();
                    else
                        bval = bdis.readByte();

                    keychars[index] = (char) bval;
                }

                String key = new String(keychars).trim();

                int chromID;
                int chromSize;
                long childOffset;

                if (isLeaf) {
                    if (isLowToHigh) {
                        chromID = lbdis.readInt();
                        chromSize = lbdis.readInt();
                    } else {
                        chromID = bdis.readInt();
                        chromSize = bdis.readInt();
                    }

                    idChromMap.put(chromID, key);
                    chromIdMap.put(key, chromID);

                    // insert leaf items
                    BPTreeLeafNodeItem leafItem = new BPTreeLeafNodeItem(++leafCount, key, chromID, chromSize);
                    thisNode.insertItem(leafItem);
                } else {
                    // get the child node pointed to in the node item
                    if (isLowToHigh)
                        childOffset = lbdis.readLong();
                    else
                        childOffset = bdis.readLong();

                    childNode = readBPTreeNode(this.fis, childOffset, thisNode, isLowToHigh);

                    // insert child node item
                    BPTreeChildNodeItem childItem = new BPTreeChildNodeItem(item, key, childNode);
                    thisNode.insertItem(childItem);
                }

                fileOffset += itemSize;
            }

        } catch (IOException ex) {
            log.error("Error reading B+ tree node " + ex);
            throw new RuntimeException("Error reading B+ tree node \n ", ex);
        }

        // success: return node
        return thisNode;
    }


}
