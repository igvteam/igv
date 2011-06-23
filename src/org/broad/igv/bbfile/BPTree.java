/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.bbfile;

import org.apache.log4j.Logger;
import org.broad.tribble.util.SeekableStream;
import org.broad.tribble.util.LittleEndianInputStream;

import java.io.DataInputStream;
import java.io.IOException;
import java.io.ByteArrayInputStream;
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

        // read in B+ tree header - verify the B+ tree info exits
        treeHeader = new BPTreeHeader(this.fis, treeOffset, isLowToHigh);

        // log error if header not found and throw exception
        if(!treeHeader.isHeaderOK()){
            int badMagic = treeHeader.getMagic();
            log.error("Error reading B+ tree header: bad magic = " + badMagic);
            throw new RuntimeException("Error reading B+ tree header: bad magic = "
                    +  badMagic);
        }

        // assign B+ tree specifications from the header
        blockSize = treeHeader.getBlockSize();
        keySize =  treeHeader.getKeySize();
        valueSize = treeHeader.getValSize();
        itemCount = treeHeader.getItemCount();

        // populate the tree - read in the nodes
        long nodeOffset = treeOffset + treeHeader.BPTREE_HEADER_SIZE;
        BPTreeNode parentNode = null;  // parent node of the root is itself, or null

        // get the root node - which recursively populates the remaining nodes
        rootNode =  readBPTreeNode(this.fis, nodeOffset, parentNode, isLowToHigh);

    }

    /*
    *   Method returns the file input stream handle
    * */
    public SeekableStream getFis() {
        return fis;
    }

    /*
    *   Method returns the B+ tree file location
    * */
    public long getBPTreeOffset() {
        return treeOffset;
    }

    /*
    *   Method returns the B+ tree header (Table E).
    * */
    public BPTreeHeader getTreeHeader(){
        return treeHeader;
    }

    /*
    *   Method returns the node block size (B+ order).
    * */
    public int getBlockSize() {
        return blockSize;
    }

    /*
    *   Method returns the chromosome name key size, which is
    *   the number of valid characters for chromosome name.
    * */
    public int getKeySize() {
        return keySize;
    }

    /*
    *   Method returns the indexing value size (currently 8).
    * */
    public int getValueSize() {
          return valueSize;
    }

    /*
    *   Method returns the number of chromosome/contig names.
    * */
    public long getItemCount() {
        return itemCount;
    }

    /*
    *   Method returns the number of nodes in the B+ tree.
    * */
    public long getNodeCount() {
        return nodeCount;
    }

    /*
    *   Method returns the root node, from which all other nodes
    *   can be extracted.
    *
    *   Returns:
    *       Root node
    * */
    public BPTreeNode getRootNode() {
        return rootNode;
    }


    Map<String, String> chromosomeKeyCache = new HashMap();
    /*
    *   Returns a search key for the mChromosome region  which  can
    *   be used to search for a corresponding section in the B+ tree.
    *
    *   According the the spec the key is the "first keySize characters of chromosome name, padded with zeroes if needed.
    * */
    public String getChromosomeKey(String chromosome) {

        String key = chromosomeKeyCache.get(chromosome);
        if(key == null) {
            char [] keyChars = new char[keySize];
            char [] chrChars = chromosome.toCharArray();
            System.arraycopy(chrChars, 0, keyChars, 0, Math.min(keySize, chrChars.length));
            key = new String(keyChars);
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
         int chromosomeID;

        // Search the B+ tree to extract the Chromosome ID.
        BPTreeNode thisNode = rootNode;

        chromosomeID = findChromosomeID(thisNode, chromKey);

        return chromosomeID;
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
         String chromKey;

        // Search the B+ tree to extract the Chromosome ID.
        BPTreeNode thisNode = rootNode;

        chromKey = findChromosomeName(thisNode, chromID);

        return chromKey;
    }

    /*
    *   Method returns all chromosome key names in B+ tree.
    *
    *   Returns:
    *   Collection of all (chromosome ID, chromosome name)entries
    * */
    public ArrayList<String> getChromosomeNames(){

        // Search the B+ tree to extract the chromosome ID.
        BPTreeNode thisNode = rootNode;

        ArrayList<String> chromosomeList = new ArrayList<String>();

        findAllChromosomeNames(thisNode, chromosomeList);

        return chromosomeList;
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
    public HashMap<Integer, String> getChromosomeIDMap(int startChromID, int endChromID){

        // Search the B+ tree to extract the chromosome ID.
        BPTreeNode thisNode = rootNode;

        HashMap<Integer, String> chromosomeIDMap = new HashMap<Integer, String>();

        findChromosomeMap(thisNode, startChromID, endChromID, chromosomeIDMap);

        return chromosomeIDMap;
    }

    // prints out the B+ Tree  nodes and leaves
    public void print() {

       // check if read in
       if(!treeHeader.isHeaderOK()){
            int badMagic = treeHeader.getMagic();
            log.error("Error reading B+ tree header: bad magic = " + badMagic);
           return;
       }

        // print B+ tree header
        treeHeader.print();

        // print  B+ tree node and leaf items - recursively
        if(rootNode != null)
            rootNode.printItems();
   }

    /*
    *   Method finds and returns the chromosome ID for the specified chromosome key.
    *
    *   Note: This method recursively calls itself, traversing the full B+ tree until
    *       either the chromosome name key is found and returns a valid chromosome ID,
    *       or exits with a -1 value.
    *
    *   Parameters:
    *       thisNode - tree node to start search
    *       chromKey - chromosome name key of valid key size.
    *
    *   Returns:
    *       Valid chromosome ID if >= 0; else -1 for not found.
    * */
    private int findChromosomeID( BPTreeNode thisNode, String chromKey){
        int chromID = -1;    // until found

        // search down the tree recursively starting with the root node
        if(thisNode.isLeaf())
        {
           int nLeaves = thisNode.getItemCount();
           for(int index = 0; index < nLeaves; ++index){
               BPTreeLeafNodeItem leaf = (BPTreeLeafNodeItem)thisNode.getItem(index);
               if(leaf == null){
                    log.error("Error finding B+ tree leaf nodes, corruption suspected");
                    throw new RuntimeException("Error reading B+ tree leaf nodes, corruption suspected");
               }

               // test chromosome key match
               if(leaf.chromKeysMatch(chromKey)){
                   chromID = leaf.getChromID();
                   break;
               }
               // else check next leaf
           }
        }
        else {
           // check all child nodes
           int nNodes = thisNode.getItemCount();
           for(int index = 0; index < nNodes; ++index){

               BPTreeChildNodeItem childItem = (BPTreeChildNodeItem)thisNode.getItem(index);              
               BPTreeNode childNode =  childItem.getChildNode();

               // check if key is in the node range
               String lowestKey = childNode.getLowestChromKey();
               String highestKey = childNode.getHighestChromKey();

               // test name key against key range
               if(chromKey.compareTo(lowestKey) >= 0
                       && chromKey.compareTo(highestKey) <= 0) {

                    // keep going until leaf items are checked
                    chromID = findChromosomeID(childNode, chromKey);

                    // check for chromKey match
                    if(chromID >= 0)
                        break;
               }
           }
        }

        return chromID;
    }

    /*
    *   Method finds and returns the chromosome name for the specified chromosome ID.
    *
    *   Parameters:
    *       thisNode - tree node to start search
    *       chromID - B+ tree chromosome ID supplied for the chromosome key
    *
    *   Returns:
    *       chromosome name if found; else a null string.
    * */
    private String findChromosomeName( BPTreeNode thisNode, int chromID){

        String chromKey = null; // mark unfound condition as an empty string

        // search down the tree recursively starting with the root node
        if(thisNode.isLeaf())
        {
           int nLeaves = thisNode.getItemCount();
           for(int index = 0; index < nLeaves; ++index){
               BPTreeLeafNodeItem leaf = (BPTreeLeafNodeItem)thisNode.getItem(index);

               if(leaf.getChromID() == chromID){ // mChromosome key match
                   chromKey = leaf.getChromKey();
                   break;
               }
               // else check next leaf
           }
        }
        else {
           // check all child nodes
           int nNodes = thisNode.getItemCount();
           for(int index = 0; index < nNodes; ++index){

               BPTreeChildNodeItem childItem = (BPTreeChildNodeItem)thisNode.getItem(index);
               BPTreeNode childNode =  childItem.getChildNode();

               // check if key is in the node range
               int lowestID = childNode.getLowestChromID();
               int highestID = childNode.getHighestChromID();

               // test chromosome ID against node ID range
               if(chromID >= lowestID && chromID <= highestID) {

                    // keep going until leaf items are checked
                    chromKey = findChromosomeName(childNode, chromID);

                    // check for chromosome ID match
                    if(chromKey != null)
                        break;
               }
           }
        }

        return chromKey;
    }

    /*
    *   Method finds and returns all chromosome names in the B+ tree.
    *
    *   Note: This method calls itself recursively until the full B+ tree is traversed.
    *
    *   Parameters:
    *       thisNode - tree node to start search
    *       chromosomeList - list of all chromosome names found.
    *
    *   Returns:
    *       Chromosome names found are added to the chromosome list passed in.
    * */
    public void findAllChromosomeNames( BPTreeNode thisNode, ArrayList<String> chromosomeList){

        // search down the tree recursively starting with the root node
        if(thisNode.isLeaf())
        {
           // add all leaf names
           int nLeaves = thisNode.getItemCount();
           for(int index = 0; index < nLeaves; ++index){

               BPTreeLeafNodeItem leaf = (BPTreeLeafNodeItem)thisNode.getItem(index);
               chromosomeList.add(leaf.getChromKey());
           }
        }
        else {
           // get all child nodes
           int nNodes = thisNode.getItemCount();
           for(int index = 0; index < nNodes; ++index){

               BPTreeChildNodeItem childItem = (BPTreeChildNodeItem)thisNode.getItem(index);
               BPTreeNode childNode = childItem.getChildNode();

               // keep going until leaf items are extracted
               findAllChromosomeNames(childNode, chromosomeList);
           }
        }
    }

    /*
    *   Method finds and returns (chromosome ID, chromosome key name) pairs for the specified ID range.
    *
    *   Parameters:
    *       thisNode - tree node to start search
    *       startChromID - starting chromosome ID for the chromosome range
    *       endChromID - ending chromosome ID for the chromosome range
    *
    *   Returns:
    *       (chromosome ID, chromosome key name) items are added to the collection passed in.
    * */
    private void findChromosomeMap( BPTreeNode thisNode, int startChromID, int endChromID,
                                        HashMap<Integer, String> chromosomeMap){
        int chromID;
        int lowestID;
        int highestID;

        // check if node is disjoint
        lowestID = thisNode.getLowestChromID();
        if(lowestID > endChromID)
            return;

        highestID = thisNode.getHighestChromID();
        if(highestID < startChromID)
            return; 

        // search down the tree recursively starting with the root node
        if(thisNode.isLeaf())
        {
           int nLeaves = thisNode.getItemCount();
           for(int index = 0; index < nLeaves; ++index){

               BPTreeLeafNodeItem leaf = (BPTreeLeafNodeItem)thisNode.getItem(index);
               chromID = leaf.getChromID();

               // check for chromosome range match
               if( chromID >= startChromID && chromID <= endChromID ){
                   chromosomeMap.put(chromID, leaf.getChromKey());
               }
               // leaf ID's are in ascending order; check for going out of range
               else if(chromID > endChromID)
                   break;
           }
        }
        else {
           // check all child nodes
           int nNodes = thisNode.getItemCount();
           for(int index = 0; index < nNodes; ++index){

               BPTreeChildNodeItem childItem = (BPTreeChildNodeItem)thisNode.getItem(index);
               BPTreeNode childNode =  childItem.getChildNode();

               // check if keys are in the node range
               lowestID = childNode.getLowestChromID();
               highestID = childNode.getHighestChromID();

               // test for chromosome range intersections
               if( lowestID <= endChromID && highestID >= startChromID )
                    findChromosomeMap(childNode, startChromID, endChromID, chromosomeMap);

               // test node ID range which is always in ascending order going out of range
               else if(lowestID > endChromID)
                   break;   //
           }
        }
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
                                      BPTreeNode parent, boolean isLowToHigh){

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

           if(isLowToHigh)
                lbdis = new LittleEndianInputStream(new ByteArrayInputStream(buffer));
           else
                bdis = new DataInputStream(new ByteArrayInputStream(buffer));

           // find node type
           if(isLowToHigh)
                type = lbdis.readByte();
           else
                type = bdis.readByte();

           // create the B+ tree node
           if(type == 1) {
               isLeaf = true;
               thisNode = new BPTreeLeafNode(++nodeCount);
           }
           else {
               isLeaf = false;
               thisNode = new BPTreeChildNode(++nodeCount);
           }

           if(isLowToHigh) {
                bval = lbdis.readByte();      // reserved - not currently used
                itemCount = lbdis.readShort();
           }
           else {
                bval = bdis.readByte();      // reserved - not currently used
                itemCount = bdis.readShort();
           }

            // Note: B+ tree node item size is the same for leaf and child items
            itemSize =  BPTREE_NODE_ITEM_SIZE + this.keySize;
            int totalSize = itemSize * itemCount;
            byte[] itemBuffer = new byte[totalSize];
            fis.readFully(itemBuffer);

            if(isLowToHigh)
                 lbdis = new LittleEndianInputStream(new ByteArrayInputStream(itemBuffer));
             else
                 bdis = new DataInputStream(new ByteArrayInputStream(itemBuffer));

            // get the node items - leaves or child nodes
            for(int item = 0; item < itemCount; ++item) {

               // always extract the key from the node format
               char[] keychars = new char[keySize];  // + 1 for 0 byte
               int index;
               for(index = 0; index < keySize; ++index) {

                    if(isLowToHigh)
                        bval = lbdis.readByte();
                    else
                        bval = bdis.readByte();

                    keychars[index] = (char)bval;
               }

               String key = new String(keychars);
                
               int chromID;
               int chromSize;
               long childOffset;

               if(isLeaf) {
                    if(isLowToHigh) {
                        chromID = lbdis.readInt();
                        chromSize = lbdis.readInt();
                    }
                    else {
                        chromID = bdis.readInt();
                        chromSize = bdis.readInt();
                    }

                    // insert leaf items
                    BPTreeLeafNodeItem leafItem = new BPTreeLeafNodeItem(++leafCount, key, chromID, chromSize);
                    thisNode.insertItem(leafItem);
               }
               else {
                   // get the child node pointed to in the node item
                   if(isLowToHigh)
                        childOffset =  lbdis.readLong();
                   else
                        childOffset =  bdis.readLong();

                   childNode = readBPTreeNode(this.fis, childOffset, thisNode, isLowToHigh);

                   // insert child node item 
                   BPTreeChildNodeItem childItem = new BPTreeChildNodeItem(item, key, childNode);
                   thisNode.insertItem(childItem);
                }

                 fileOffset += itemSize;
           }

        }catch(IOException ex) {
           log.error("Error reading B+ tree node " + ex);
           throw new RuntimeException("Error reading B+ tree node \n ", ex);
        }

        // success: return node
        return thisNode;
   }


}