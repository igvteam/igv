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

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: martind
 * Date: Dec 17, 2009
 * Time: 12:28:30 PM
 * To change this template use File | Settings | File Templates.
 */
/*
*   Container class for B+ Tree Leaf Node.
*
*   Note: Key is a property of node items and BPTreeNode methods
*       getLowestKeyItem() and getHighestKeyItem() can be used
*       to check node key range.
* */
public class BPTreeLeafNode implements BPTreeNode{

    private static Logger log = Logger.getLogger(BPTreeLeafNode.class);
    private final boolean isLeafNode = true;

    private long nodeIndex;    // index for node in B+ tree organization
    String lowestChromKey;     // lowest chromosome/contig key name
    String highestChromKey;    // highest chromosome/contig key name
    int lowestChromID;         // lowest chromosome ID corresponds to lowest key
    int highestChromID;        // highest chromosome ID corresponds to highest key
    private ArrayList<BPTreeLeafNodeItem> leafItems; // array for leaf items

    /*
    *   Constructor for the B+ tree leaf (terminal) node.
    *
    *   Parameters:
    *       nodeIndex - index assigned to the node
    *       parent - parent node (object)
    *
    *   Note: Inserted leaf items contain associated name key/chromosome ID.
    * */
    public BPTreeLeafNode(long nodeIndex){

        this.nodeIndex = nodeIndex;
        leafItems = new ArrayList<BPTreeLeafNodeItem>();
    }

    /*
    *   Method returns the node index in the B+ tree organization.
    *
    *   Returns:
    *       node index in B+ tree
    * */
    public long getNodeIndex(){
         return nodeIndex;
    }

    /*
    *   Method identifies the node as a leaf node or a child (non-leaf) node.
    *
    *   Returns:
    *       true, if leaf node; false if child node
    * */
    public boolean isLeaf() {
        return isLeafNode;
    }

    /*
    *   Method inserts the node item appropriate to the item's key value.
    *
    *   Returns:
    *       Node item inserted successfully.
    * */
    public boolean insertItem(BPTreeNodeItem item){

         // Quick implementation: assumes all keys are inserted in rank order
        // todo: verify if need to compare key and insert at rank location
        leafItems.add((BPTreeLeafNodeItem)item );

        // Note: assumes rank order insertions
        if(leafItems.size() == 1 ){
            lowestChromKey = item.getChromKey();
            lowestChromID = ((BPTreeLeafNodeItem)item).getChromID();
        }
        else {
           highestChromKey = item.getChromKey();
           highestChromID = ((BPTreeLeafNodeItem)item).getChromID();
        }

        // success
        return true;
    }

    /*
    *   Method deletes the node item appropriate to the item's index.
    *
    *   Returns:
    *       Node item deleted successfully.
    * */
    public boolean deleteItem(int index){

        // unacceptable index
        if(index < 0 || index >= getItemCount())
           return false;

        leafItems.remove(index);
        return true;  // success
    }

    /*
    *   Method returns the number of items assigned to the node.
    *
    *   Returns:
    *       Count of node items contained in the node
    * */
    public int getItemCount() {
        return leafItems.size();
    }

    /*
    *   Method returns the indexed node item.
    *
    *   Returns:
    *       Indexed node item.
    * */
    public  BPTreeNodeItem getItem(int index){
        if(getItemCount() > 0 && index < getItemCount())
            return leafItems.get(index);
        else
            return null;
    }

    /*
    *   Method returns the lowest key value belonging to the node.
    *
    *   Returns:
    *       Lowest key contig/chromosome name value
    * */
    public  String getLowestChromKey(){
       if(leafItems.size() > 0)
           return lowestChromKey;
       else
           return null;
    }

    /*
    *   Method returns the highest key value belonging to the node.
    *
    *   Returns:
    *       Highest key contig/chromosome name value
    * */
    public  String getHighestChromKey(){
       if(leafItems.size() > 0)
           return highestChromKey;
       else
           return null;
    }

     /*
    *   Method returns the lowest chromosome ID belonging to the node.
    *
    *   Returns:
    *       Lowest key contig/chromosome ID; or -1 if no node items
    * */
    public  int getLowestChromID(){
        if(leafItems.size() > 0)
            return lowestChromID;
        else
            return -1;
    }

    /*
    *   Method returns the highest chromosome ID belonging to the node.
    *
    *   Returns:
    *       Highest key contig/chromosome ID; or -1 if no node items
    * */
    public  int getHighestChromID(){
        if(leafItems.size() > 0)
            return highestChromID;
        else
            return -1;
    }

    /*
    *   Method prints the nodes items and sub-node items.
    *       Node item deleted successfully.
    * */
    public void printItems(){
        int  itemCount = getItemCount();

        log.debug("Leaf node " + nodeIndex +  "contains " + itemCount + " leaf items:");
        for(int item = 0; item < itemCount; ++item){
            leafItems.get(item).print();
        }
    }

    // ************** BPTreeLeafNode specific methods ***********
    /*
    *   Method returns all leaf items mContained by this leaf node.
    *
    *   Returns:
    *       List of leaf items contained by this node
    * */
    public  ArrayList<BPTreeLeafNodeItem> getLeafItems() {
        return leafItems;
    }

}