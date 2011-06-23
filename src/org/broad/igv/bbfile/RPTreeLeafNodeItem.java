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
 * Time: 11:21:49 PM
 * To change this template use File | Settings | File Templates.
 */
/*
    Container class for R+ tree leaf node data locator.
*
*   Note: Determination of data item as  BigWig data or BigBed data
*           depends on whether the file is BigWig of Table J format
*           or BigBed of Tble I format.
 */
public class RPTreeLeafNodeItem extends RPChromosomeRegion implements RPTreeNodeItem {

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
                              int endChromID, int endBase, long dataOffset, long dataSize) {
        super(startChromID, startBase, endChromID, endBase);
        this.dataOffset = dataOffset;
        this.dataSize = dataSize;
    }

 
    public RPChromosomeRegion getChromosomeBounds() {
        return this;
    }

    public void print() {

        log.debug("R+ tree leaf node data item ");
        log.debug("StartChromID = " + getStartChromID());
        log.debug("StartBase = " + getStartBase());
        log.debug("EndChromID = " + getEndChromID());
        log.debug("EndBase = " + getEndBase());

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
