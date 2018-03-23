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
public class RPTreeChildNodeItem extends RPTreeNodeItem {

    private static Logger log = Logger.getLogger(RPTreeChildNodeItem.class);

    private RPTreeNode childNode;  // child node assigned to node item
    private RPTreeNodeProxy childNodeProxy;

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

        super(new RPChromosomeRegion(startChromID, startBase, endChromID, endBase));
        this.childNode = childNode;
    }


    public RPTreeChildNodeItem(int startChromID, int startBase,
                               int endChromID, int endBase, RPTreeNodeProxy childNodeProxy) {
        super(new RPChromosomeRegion(startChromID, startBase, endChromID, endBase));
        this.childNodeProxy = childNodeProxy;
    }

    public RPTreeNode getChildNode() {

        if (childNode == null) {
            RPTreeNodeProxy proxy = childNodeProxy;
            childNode = RPTree.readRPTreeNode(proxy.fis, proxy.fileOffset, proxy.isLowToHigh, true);
        }

        return childNode;
    }

    public void print() {

        super.print();

        // child node specific entries
        childNode.printItems();
    }

}

