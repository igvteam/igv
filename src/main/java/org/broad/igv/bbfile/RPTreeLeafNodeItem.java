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
 * @author Martin Decautis
 */

/**
 * Container class for R+ tree leaf node data locator.
 * <p/>
 * Note: Determination of data item as  BigWig data or BigBed data
 * depends on whether the file is BigWig of Table J format
 * or BigBed of Tble I format.
 */
public class RPTreeLeafNodeItem extends RPTreeNodeItem {

    private static Logger log = Logger.getLogger(RPTreeLeafNodeItem.class);

    private long dataOffset;      // file offset to data item
    private long dataSize;        // size of data item

    /*  Constructor for leaf node items.
    *
    *   Parameters:
    *       itemIndex - index of item belonging to a leaf node
    *       startChromID - starting chromosome/contig for item
    *       startBase - starting base for item
    *       endChromID - ending chromosome/contig for item
    *       endBase - ending base for item
    *       dataOffset - file location for leaf chromosome/contig data
    *       dataSize - size of (compressed) leaf data region in bytes
    *
    * */

    public RPTreeLeafNodeItem(int startChromID, int startBase,
                              int endChromID, int endBase,
                              long dataOffset, long dataSize) {
        super(new RPChromosomeRegion(startChromID, startBase, endChromID, endBase));
        this.dataOffset = dataOffset;
        this.dataSize = dataSize;
    }

    public void print() {

        log.debug("R+ tree leaf node data item ");
        super.print();

        // leaf node specific entries
        log.debug("DataOffset = " + dataOffset);
        log.debug("DataSize = " + dataSize);
    }

    // *** RPTreeLeafNodeItem specific methods ***

    public long getDataOffset() {
        return dataOffset;
    }

    public long geDataSize() {
        return dataSize;
    }

}
