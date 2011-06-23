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

/**
 * Created by IntelliJ IDEA.
 * User: martind
 * Date: Dec 20, 2009
 * Time: 10:50:26 PM
 * To change this template use File | Settings | File Templates.
 */
/*
    Container class for B+ Tree child node format
 */
public class BPTreeChildNodeItem implements BPTreeNodeItem {

    private static Logger log = Logger.getLogger(BPTreeChildNodeItem.class);
    private final boolean isLeafItem = false;
    private long itemIndex;     // item index in child node list

    // B+ Tree Child Node Item entities - BBFile Table H
    // Note the childOffset entity is replaced with the actual instance
    // of the child note it points to in the file.
    private String chromKey;   // mChromosome/contig name; of keysize chars
    private BPTreeNode childNode;  // child node

    /*
    *   Constructs a B+ tree child node item with the supplied information.
    *
    *   Parameters:
    *       itemIndex - node item index
    *       chromKey - chromosome name key
    *       childNode - assigned child node object
    * */
    public BPTreeChildNodeItem(int itemIndex, String chromKey, BPTreeNode childNode){
        this.itemIndex =  itemIndex;
        this.chromKey = chromKey;
        this.childNode = childNode;
    }

    /*
    *   Method returns the index assigned to this node item.
    *
    *   Returns:
    *       index assigned to this node item
    * */
     public long getItemIndex() {
        return itemIndex;
    }

    /*
    *   Method returns if this node is a leaf item.
    *
    *   Returns:
    *       false because node is a child (non-leaf) item
    * */
     public boolean isLeafItem() {
        return isLeafItem;
    }

    /*
    *   Method returns the chromosome name key  assigned to this node item.
    *
    *   Returns:
    *       chromosome name key assigned to this node item
    * */
     public String getChromKey() {
        return chromKey;
    }

    /*
    *   Method compares supplied chromosome key with leaf node key.
    *
    *   Parameters:
    *       chromKey - chromosome name ley to compare
    *
    *   Returns:
    *       true, if keys are equal; false if keys are different
    * */
    public boolean chromKeysMatch(String chromKey) {
        String thisKey = this.chromKey;
        String thatKey = chromKey;

        // Note: must have the same length to compare chromosome names
        int thisKeyLength = thisKey.length();
        int thatKeyLength = thatKey.length();

        // check if need to truncate the larger string
        if(thisKeyLength > thatKeyLength)
            thisKey = thisKey.substring(0,thatKeyLength);
        else if(thatKeyLength > thisKeyLength)
            thatKey = thatKey.substring(0,thisKeyLength);

        if (thisKey.compareTo(thatKey) == 0)
            return true;
        else
            return false;
    }

    public void print() {

        log.debug("B+ Tree child node " + itemIndex);
        log.debug("Key value = " + chromKey);

        // recursively print chid node items
        childNode.printItems();
   }

    // BPTreeLeafNodeItem specific methods
    public BPTreeNode getChildNode() {
        return childNode;
    }

}
