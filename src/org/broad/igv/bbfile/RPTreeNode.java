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

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: martind
 * Date: Jan 14, 2010
 * Time: 11:16:03 AM
 * To change this template use File | Settings | File Templates.
 */
public class RPTreeNode {

    boolean leaf;
    protected RPChromosomeRegion chromosomeBounds;  // chromosome bounds for entire node
    protected ArrayList<RPTreeNodeItem> items; // array for child items

    public RPTreeNode(boolean leaf) {
        this.leaf = leaf;
        items = new ArrayList<RPTreeNodeItem>();
    }

    // Identifies the node as a leaf node or a child (non-leaf) node.
    public boolean isLeaf() {
        return leaf;
    }

    public RPChromosomeRegion getChromosomeBounds() {
        return chromosomeBounds;
    }

    public int compareRegions(RPChromosomeRegion chromosomeRegion) {
        return chromosomeBounds.compareRegions(chromosomeRegion);
    }

    public int getItemCount() {
        return items.size();
    }

    public RPTreeNodeItem getItem(int index) {
        if (index < 0 || index >= items.size())
            return null;
        else {
            return items.get(index);
        }
    }

    public void insertItem(RPTreeNodeItem item) {

        RPTreeNodeItem newItem = item;

        // Quick implementation: assumes all keys are inserted in rank order
        // todo: or compare key and insert at rank location
        items.add(newItem);

        // Update node bounds or start node chromosome bounds with first entry
        if (chromosomeBounds == null) {
            chromosomeBounds = new RPChromosomeRegion(newItem.getChromosomeBounds());
        } else {
            chromosomeBounds = chromosomeBounds.getExtremes(newItem.getChromosomeBounds());
        }

    }

    public void printItems() {

        for (int item = 0; item < items.size(); ++item) {
            items.get(item).print();
        }
    }
}
