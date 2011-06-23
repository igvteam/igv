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

/**
 * Created by IntelliJ IDEA.
 * User: martind
 * Date: Jan 6, 2010
 * Time: 4:05:31 PM
 * To change this template use File | Settings | File Templates.
 */
/*
*   RPTree class will construct a R+ tree from a binary Bed/Wig BBFile.
*   (or by insertion of tree nodes - TBD see insert method)
*
*   1) RPTree will first read in the R+ tree header with RPTreeHeader class.
*
*   2) Starting with the root node, the readRPTreeNode method will read in the
*   node format, determine if the node contains child nodes (isLeaf = false)
*   or leaf items (isLeaf = true).
*
*   3) If the node is a leaf node, all leaf items are read in to the node's leaf array.
*
*   4) If node is a child node, readRPTreeNode will be called recursively,
*   until a leaf node is encountered, where step 3 is performed.
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
public class RPTree {

    private static Logger log = Logger.getLogger(RPTree.class);

    public static final int RPTREE_NODE_FORMAT_SIZE = 4;       // node format size
    public static final int RPTREE_NODE_LEAF_ITEM_SIZE = 32;   // leaf item size
    public static final int RPTREE_NODE_CHILD_ITEM_SIZE = 24;  // child item size

    // R+ tree access variables   - for reading in R+ tree nodes from a file
    private int uncompressBuffSize;    // decompression buffer size; or 0 for uncompressed data
    private boolean isLowToHigh;       // binary data low to high if true; else high to low
    private long rpTreeOffset;         // file offset to the R+ tree

    // R+ tree index header - Table K
    private RPTreeHeader rpTreeHeader; // R+ tree header (Table K for BBFile)

    // R+ tree bounds
    private RPChromosomeRegion chromosomeBounds;  // R+ tree's bounding chromosome region

    // R+ tree nodal variables
    private int order;         // R+ tree order: maximum number of leaves per node
    private RPTreeNode rootNode;  // root node for R+ tree
    private long nodeCount;        // number of nodes defined in the R+ tree
    private long leafCount;        // number of leaves in the R+ tree


    /*
   * Constructor for reading in a B+ tree from a BBFile/input stream.
   * */
    /*
    *   Constructor for R+ chromosome data locator tree
    *
    *   Parameters:
    *       fis - file input stream handle
    *       fileOffset - location for R+ tree header
    *       isLowToHigh - binary values are low to high if true; else high to low
    *       uncompressBuffSize - buffer size for decompression; else 0 for uncompressed data
    * */

    public RPTree(SeekableStream fis, long fileOffset, boolean isLowToHigh, int uncompressBuffSize, boolean forceDescend) {

        // save the seekable file handle  and B+ Tree file offset
        // Note: the offset is the file position just after the B+ Tree Header
        // mBBFis = fis;
        rpTreeOffset = fileOffset;
        this.uncompressBuffSize = uncompressBuffSize;
        this.isLowToHigh = isLowToHigh;

        // read in R+ tree header - verify the R+ tree info exits
        rpTreeHeader = new RPTreeHeader(fis, rpTreeOffset, isLowToHigh);

        // log error if header not found and throw exception
        if (!rpTreeHeader.isHeaderOK()) {
            int badMagic = rpTreeHeader.getMagic();
            log.error("Error reading R+ tree header: bad magic = " + badMagic);
            throw new RuntimeException("Error reading R+ tree header: bad magic = " + badMagic);
        }

        // assigns R+ tree organization from the header
        order = rpTreeHeader.getBlockSize();
        chromosomeBounds = new RPChromosomeRegion(rpTreeHeader.getStartChromID(), rpTreeHeader.getStartBase(),
                rpTreeHeader.getEndChromID(), rpTreeHeader.getEndBase());

        // populate the tree - read in the nodes
        long nodeOffset = rpTreeOffset + rpTreeHeader.getHeaderSize();
        RPTreeNode parentNode = null;      // parent node of the root is itself, or null

        // start constructing the R+ tree - get the root node
        rootNode = readRPTreeNode(fis, nodeOffset, isLowToHigh, forceDescend);
    }

    /*
     *  Constructs an R+ Tree which conforms to the supplied information
     *      order -  the items per node factor, sometimes called the m factor,
     *      where any node must have at least m/2 items and no more than m items.
     *      keySize - the number of significant bytes in a item key.
     * */

    public RPTree(int order) {

        // R+ tree node specification
        this.order = order;

        // Note: acknowledge no bounds specified as a null  object
        chromosomeBounds = null;

    }

    /*
    *   Method returns size of buffer required for data decompression.
    *
    *   Returns:
    *       Data decompression buffer size in bytes; else 0 for uncompressed data in file.
    * */

    public int getUncompressBuffSize() {
        return uncompressBuffSize;
    }

    /*
    *   Method returns if file contains formatted data in low to high byte order.
    *
    *   Returns:
    *       Returns true if data ordered in low to high byte order; false if high to low.
    * */

    public boolean isIsLowToHigh() {
        return isLowToHigh;
    }

    /*
    *   Method returns the R+ tree order, the maximum number of leaf items per node.
    *
    *   Returns:
    *       Maximum number of leaf items per node..
    * */

    public int getOrder() {
        return order;
    }

    /*
    *   Method returns the R+ tree index header.
    *
    *   Returns:
    *       R+ tree index header.
    * */

    public RPTreeHeader getRPTreeHeader() {
        return rpTreeHeader;
    }

    /*
    *   Method returns the total number of chromosomes or contigs in the R+ tree.
    *
    *   Returns:
    *       Total number of chromosomes or contigs in the R+ tree.
    * */

    public long getItemCount() {
        return rpTreeHeader.getItemCount();
    }

    /*
    *   Method returns the chromosome bounding region for all R+ tree data.
    *
    *   Returns:
    *       chromosome bounding region for all R+ tree data
    * */

    public RPChromosomeRegion getChromosomeBounds() {
        return chromosomeBounds;
    }

    /*
    *   Method returns the total node count for the R+ tree.
    *
    *   Returns:
    *       Node count for R+ tree; or 0 if tree was constructed without nodes.
    * */

    public long getNodeCount() {
        return nodeCount;
    }

    /*
    *   Method finds the bounding chromosome region in R+ tree for a chromosome ID range.
    *
    *   Parameters:
    *       startChromID - start chromosome for the region
    *       endChromID - end chromosome for the region
    *
    *   Returns:
    *       Region which bounds the extremes of chromosome ID range
    * */

    public RPChromosomeRegion getChromosomeRegion(int startChromID, int endChromID) {

        RPChromosomeRegion region;

        // Search the R+ tree to extract the chromosome region.
        RPTreeNode thisNode = rootNode;
        RPChromosomeRegion seedRegion = null;  // null until a chromosome match

        region = findChromosomeRegion(thisNode, startChromID, endChromID, seedRegion);

        return region;
    }

    /*
    *   Method returns list of all chromosome regions found for the chromosome ID range.
    *
    *   Returns:
    *       List of all chromosome regions in the chromosome ID range.
    * */

    public ArrayList<RPChromosomeRegion> getAllChromosomeRegions() {

        // Search the R+ tree to extract the chromosome regions
        RPTreeNode thisNode = rootNode;

        ArrayList<RPChromosomeRegion> regionList = new ArrayList<RPChromosomeRegion>();

        findAllChromosomeRegions(thisNode, regionList);

        return regionList;
    }

    /*
    *   Method extracts a hit list of chromosome data file locations for a specified chromosome region.
    *
    *   Parameters:
    *       chromosomeRegion - chromosome region for feature extraction consists of:
    *           startChromID - start chromosome ID for region
    *           mStartBase - starting base for data extraction
    *           endChromID - end chromosome ID for region
    *           mEndBase - ending base for data extraction
    *       contained - if true indicates all returned data must be
    *           completely contained within the extraction region;
    *           else if false, returns all intersecting region features
    *
    *   Note: The selection region will be limited to accommodate  maxLeafHits; which terminates
    *       selection at the leaf node at which maxLeafHits is reached. Total number of selected
    *       items may exceed maxLeafHits, but only by the number of leaves in the cutoff leaf node.
    *
    *   Returns:
    *       List of chromosome leaf items which identify file locations for bed data
    *       of a chromosome region, or a sub-region subject to maxLeafHits.
    *
    *       Check returned leaf item bounds for cutoff limits on selection region due to maxLeafHits.
    * */

    public ArrayList<RPTreeLeafNodeItem> getChromosomeDataHits(RPChromosomeRegion selectionRegion,
                                                               boolean contained) {

        ArrayList<RPTreeLeafNodeItem> leafHitItems = new ArrayList<RPTreeLeafNodeItem>();

        // check for valid selection region - return empty collection if null
        if (selectionRegion == null)
            return leafHitItems;

        // limit the hit list size
        /*
        if(maxLeafHits > 0)
            mMaxLeafHits = maxLeafHits;
        else
            mMaxLeafHits = mRPTreeHeader.getBlockSize();
        */


        findChromosomeRegionItems(rootNode, selectionRegion, leafHitItems);

        return leafHitItems;
    }

    // prints out the R+ tree  header, nodes, and leaves

    public void print() {

        // check if read in
        if (!rpTreeHeader.isHeaderOK()) {
            int badMagic = rpTreeHeader.getMagic();
            log.error("Error reading R+ tree header: bad magic = " + badMagic);
            return;
        }

        // print R+ tree header
        rpTreeHeader.print();

        // print  R+ tree node and leaf items - recursively
        if (rootNode != null)
            rootNode.printItems();

    }

    /*
    *   Method finds and returns the bounding chromosome region for the specified
    *   chromosome ID range.
    *
    *   Parameters:
    *       thisNode - tree node to start search
    *       startChromID - start chromosome ID for region
    *       endChromID - end chromosome ID for region
    *       region  - leaf region contains extremes for given chromosome ID range
    *
    *   Note: region grows recursively to match the extremes found for the
    *   specified chromosome ID range.  Starting base comes from the startChromID
    *   match and ending base comes form the endChromID match.
    *
    *   Returns:
    *       Chromosome region if found in the R+ tree node passed in;
    *       else null region
    * */

    private RPChromosomeRegion findChromosomeRegion(RPTreeNode thisNode,
                                                    int startChromID, int endChromID, RPChromosomeRegion region) {

        int hitValue;
        RPChromosomeRegion bounds;

        // search down the tree recursively starting with the root node
        if (thisNode.isLeaf()) {
            int nLeaves = thisNode.getItemCount();
            for (int index = 0; index < nLeaves; ++index) {

                RPTreeLeafNodeItem leaf = (RPTreeLeafNodeItem) thisNode.getItem(index);

                // get leaf region bounds
                bounds = leaf.getChromosomeBounds();

                // test this leaf's chromosome ID's for chromosome hit, then include its base bounds
                if (startChromID >= bounds.getStartChromID() && startChromID <= bounds.getEndChromID() ||
                        endChromID >= bounds.getStartChromID() && endChromID <= bounds.getEndChromID()) {

                    // Note: need a start region before comparing other regions for extremes
                    if (region == null)
                        region = new RPChromosomeRegion(bounds); // seed extreme region
                    else
                        region = region.getExtremes(bounds); // update seed extreme region
                }
            }
        } else {
            // check all child nodes
            int nNodes = thisNode.getItemCount();
            for (int index = 0; index < nNodes; ++index) {

                RPTreeChildNodeItem childItem = (RPTreeChildNodeItem) thisNode.getItem(index);

                // get bounding region and compare chromosome ID's
                bounds = childItem.getChromosomeBounds();

                // test node chromosome ID range for any leaf hits for either startChromID or endChromID
                if (startChromID >= bounds.getStartChromID() && startChromID <= bounds.getEndChromID() ||
                        endChromID >= bounds.getStartChromID() && endChromID <= bounds.getEndChromID()) {

                    RPTreeNode childNode = childItem.getChildNode();
                    region = findChromosomeRegion(childNode, startChromID, endChromID, region);
                }

                // check next node
            }
        }

        return region;
    }

    /*
    *   Method finds and returns all chromosome regions in the R+ chromosome data tree.
    *
    *   parameters:
    *       thisNode - tree node to start search
    *       chromosomeList - list of all chromosome names found.
    *
    *   Returns:
    *       Adds chromosome regions if found in the chromosome region list passed in.
    * */

    private void findAllChromosomeRegions(RPTreeNode thisNode,
                                          ArrayList<RPChromosomeRegion> regionList) {

        // search down the tree recursively starting with the root node
        if (thisNode.isLeaf()) {
            int nLeaves = thisNode.getItemCount();
            for (int index = 0; index < nLeaves; ++index) {

                RPTreeLeafNodeItem leaf = (RPTreeLeafNodeItem) thisNode.getItem(index);

                // add all leaf regions
                RPChromosomeRegion region = leaf.getChromosomeBounds();
                regionList.add(region);
            }
        } else {
            // get all child nodes
            int nNodes = thisNode.getItemCount();
            for (int index = 0; index < nNodes; ++index) {

                RPTreeChildNodeItem childItem = (RPTreeChildNodeItem) thisNode.getItem(index);
                RPTreeNode childNode = childItem.getChildNode();

                findAllChromosomeRegions(childNode, regionList);
            }
        }

    }

    /*
    *   Method returns an array of chromosome leaf items for the chromosome test region.
    *
    *   Note: At the leaf item level, any hit is valid. Use of contained is exercised
    *       when the leaf item data region is examined.
    *   Parameters:
    *
    *       thisNode - tree node to start search
    *       testRegion - bounding region for feature extraction consists of:
    *           start chromosome ID for region
    *           starting base for data extraction
    *           end chromosome ID for region
    *           ending base for data extraction
    *       leafHitItems - array containing previous leaf hit items
    *
    *   Note: leaf hit items will be limited to leaves in the current leaf node,
    *       once the maximum number of leaf hits mMaxLeafHits is reached.
    * 
    *   Returns:
    *       ArrayList of leaf hit items containing updated hit regions and data offsets;
    *       else an empty array if hit regions not found.
    * */

    private void findChromosomeRegionItems(RPTreeNode thisNode, RPChromosomeRegion selectionRegion,
                                           ArrayList<RPTreeLeafNodeItem> leafHitItems) {

        int hitValue;

        // check for valid selection region - ignore request if null
        if (selectionRegion == null)
            return;

        // check if node is disjoint
        hitValue = thisNode.compareRegions(selectionRegion);
        if (Math.abs(hitValue) >= 2)
            return;

        // search down the tree recursively starting with the root node
        if (thisNode.isLeaf()) {
            int nLeaves = thisNode.getItemCount();
            for (int index = 0; index < nLeaves; ++index) {
                RPTreeLeafNodeItem leafItem = (RPTreeLeafNodeItem) thisNode.getItem(index);

                // compute the region hit value
                hitValue = leafItem.compareRegions(selectionRegion);

                // select contained or intersected leaf regions - item selection is by iterator
                if (Math.abs(hitValue) < 2) {
                    leafHitItems.add(leafItem);
                }

                // ascending regions will continue to be disjoint so terminate nodal search
                else if (hitValue > 1)
                    break;

                // check next leaf
            }
        } else {
            // check all child nodes
            int nNodes = thisNode.getItemCount();
            for (int index = 0; index < nNodes; ++index) {
                RPTreeChildNodeItem childItem = (RPTreeChildNodeItem) thisNode.getItem(index);

                // check for region intersection at the node level
                hitValue = childItem.compareRegions(selectionRegion);

                // test this node and get any leaf hits; intersections and containing
                if (Math.abs(hitValue) < 2) {
                    RPTreeNode childNode = childItem.getChildNode();
                    findChromosomeRegionItems(childNode, selectionRegion, leafHitItems);
                }

                // ascending regions will continue to be disjoint so terminate nodal search
                else if (hitValue > 1)
                    break;
            }
        }
    }

    /*
    *   Method reads in the R+ tree nodes recursively.
    *
    *   Note: If node is a child node, the node is examined recursively,
    *       until the leaves are found.
    *
    *   Parameters:
    *       fis - file input stream handle
    *       fileOffset - file location for node specification (Table L)
    *       parent - parent node of this node
    *       isLowToHigh - indicates formatted data is low to high byte order if true;
    *           else is high to low byte order
    *
    *   Returns:
    *       A tree node, for success, or null for failure to find the node information.

    * */

    static RPTreeNode readRPTreeNode(SeekableStream fis, long fileOffset, boolean isLowToHigh, boolean forceDescend) {

        LittleEndianInputStream lbdis = null; // low o high byte stream reader
        DataInputStream bdis = null;    // high to low byte stream reader

        byte[] buffer = new byte[RPTREE_NODE_FORMAT_SIZE];
        RPTreeNode thisNode = null;

        try {

            // Read node format into a buffer
            fis.seek(fileOffset);
            fis.readFully(buffer);

            if (isLowToHigh)
                lbdis = new LittleEndianInputStream(new ByteArrayInputStream(buffer));
            else
                bdis = new DataInputStream(new ByteArrayInputStream(buffer));

            // find node type
            byte type;
            if (isLowToHigh)
                type = lbdis.readByte();
            else
                type = bdis.readByte();


            boolean isLeaf;
            int itemSize;
            if (type == 1) {
                isLeaf = true;
                itemSize = RPTREE_NODE_LEAF_ITEM_SIZE;
                thisNode = new RPTreeLeafNode();
            } else {
                isLeaf = false;
                itemSize = RPTREE_NODE_CHILD_ITEM_SIZE;
                thisNode = new RPTreeChildNode();
            }
            //nodeCount++;

            int itemCount;
            if (isLowToHigh) {
                lbdis.readByte();          // reserved - not currently used
                itemCount = lbdis.readShort();
            } else {
                bdis.readByte();          // reserved - not currently used
                itemCount = bdis.readShort();
            }

            int itemBlockSize = itemCount * itemSize;
            buffer = new byte[itemBlockSize];            // allocate buffer for item sisze
            fis.readFully(buffer);
            if (isLowToHigh)
                lbdis = new LittleEndianInputStream(new ByteArrayInputStream(buffer));
            else
                bdis = new DataInputStream(new ByteArrayInputStream(buffer));

            // get the node items - leaves or child nodes
            int startChromID, endChromID;
            int startBase, endBase;
            for (int item = 0; item < itemCount; ++item) {


                // always extract the bounding rectangle
                if (isLowToHigh) {
                    startChromID = lbdis.readInt();
                    startBase = lbdis.readInt();
                    endChromID = lbdis.readInt();
                    endBase = lbdis.readInt();
                } else {
                    startChromID = bdis.readInt();
                    startBase = bdis.readInt();
                    endChromID = bdis.readInt();
                    endBase = bdis.readInt();
                }

                if (isLeaf) {
                    long dataOffset;
                    long dataSize;
                    if (isLowToHigh) {
                        dataOffset = lbdis.readLong();
                        dataSize = lbdis.readLong();
                    } else {
                        dataOffset = bdis.readLong();
                        dataSize = bdis.readLong();
                    }

                    thisNode.insertItem(new RPTreeLeafNodeItem(startChromID, startBase, endChromID, endBase,
                            dataOffset, dataSize));
                } else {
                    // get the child node pointed to in the node item
                    long nodeOffset;
                    if (isLowToHigh) {
                        nodeOffset = lbdis.readLong();
                    } else {
                        nodeOffset = bdis.readLong();
                    }

                    // Recursive call to read next child node
                    // The test on chromIds is designed to stop the descent when the tree reaches the level of an
                    // individual chromosome.  These are loaded later on demand.
                    RPTreeNode childNode;
                    if (startChromID != endChromID || forceDescend) {
                        childNode = readRPTreeNode(fis, nodeOffset, isLowToHigh, forceDescend);
                    } else {
                        childNode = new RPTreeNodeProxy(fis, nodeOffset, isLowToHigh, startChromID);
                    }


                    // insert child node item
                    thisNode.insertItem(new RPTreeChildNodeItem(startChromID, startBase, endChromID,
                            endBase, childNode));


                }

                fileOffset += itemSize;
            }

        } catch (IOException ex) {
            log.error("Error reading in R+ tree nodes: " + ex);
            throw new RuntimeException("Error reading R+ tree nodes: \n", ex);
        }

        // return success
        return thisNode;
    }
}
