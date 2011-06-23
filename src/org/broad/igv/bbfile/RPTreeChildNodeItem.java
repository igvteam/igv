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
 * Date: Jan 6, 2010
 * Time: 4:35:42 PM
 * To change this template use File | Settings | File Templates.
 */

/*
    Container class for R+ Tree Child format
 */
public class RPTreeChildNodeItem implements RPTreeNodeItem {

    private static Logger log = Logger.getLogger(RPTreeChildNodeItem.class);

    // R+ child (non-leaf) node item entries: BBFile Table N
    private RPChromosomeRegion chromosomeBounds; // chromosome bounds for item
    private RPTreeNode childNode;  // child node assigned to node item

    /*  Constructor for child node items.
    *
    *   Parameters:
    *       itemIndex - index of item belonging to a child node
    *       startChromID - starting chromosome/contig for item
    *       startBase - starting base for item
    *       endChromID - ending chromosome/contig for item
    *       endBase - ending base for item
    *       childNode - child node item assigned to child node
    *
    * */

    public RPTreeChildNodeItem(int startChromID, int startBase,
                               int endChromID, int endBase, RPTreeNode childNode) {


        chromosomeBounds = new RPChromosomeRegion(startChromID, startBase, endChromID, endBase);
        this.childNode = childNode;
    }


    public RPTreeChildNodeItem(int startChromID, int startBase,
                               int endChromID, int endBase, long childDataOffset) {


        chromosomeBounds = new RPChromosomeRegion(startChromID, startBase, endChromID, endBase);
    }


    public RPChromosomeRegion getChromosomeBounds() {
        return chromosomeBounds;
    }

    public RPTreeNode getChildNode() {

        if (childNode instanceof RPTreeNodeProxy) {
            RPTreeNodeProxy proxy = (RPTreeNodeProxy) childNode;
             childNode = RPTree.readRPTreeNode(proxy.fis, proxy.fileOffset, proxy.isLowToHigh, true);
        }

        return childNode;
    }

    public int compareRegions(RPChromosomeRegion chromosomeRegion) {

        int value = chromosomeBounds.compareRegions(chromosomeRegion);
        return value;
    }

    public void print() {

        log.debug("Child node item :\n");
        log.debug(" StartChromID = " + chromosomeBounds.getStartChromID() + "\n");
        log.debug(" StartBase = " + chromosomeBounds.getStartBase() + "\n");
        log.debug(" EndChromID = " + chromosomeBounds.getEndChromID() + "\n");
        log.debug(" EndBase = " + chromosomeBounds.getEndBase() + "\n");

        // child node specific entries
        childNode.printItems();
    }

}

