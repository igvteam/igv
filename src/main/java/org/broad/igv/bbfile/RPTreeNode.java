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
