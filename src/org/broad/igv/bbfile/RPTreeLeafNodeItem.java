/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not
 * responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL), Version 2.1 which is
 * available at http://www.opensource.org/licenses/lgpl-2.1.php.
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
