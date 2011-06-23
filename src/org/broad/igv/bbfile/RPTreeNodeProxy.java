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

import org.broad.tribble.util.SeekableStream;

/**
 * @author jrobinso
 * @date Jun 22, 2011
 */
public class RPTreeNodeProxy implements RPTreeNode {

    public SeekableStream fis;
    public long fileOffset;
    public boolean isLowToHigh;

    // For debugging
    int chromId;

    public RPTreeNodeProxy(SeekableStream fis, long fileOffset, boolean lowToHigh, int chromId) {
        this.fis = fis;
        this.fileOffset = fileOffset;
        isLowToHigh = lowToHigh;
        this.chromId = chromId;
    }

    public boolean isLeaf() {
        throw new UnsupportedOperationException("Not implemented -- this should never be called on a node proxy");
    }

    public RPChromosomeRegion getChromosomeBounds() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public int compareRegions(RPChromosomeRegion chromosomeRegion) {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public int getItemCount() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public RPTreeNodeItem getItem(int index) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public boolean insertItem(RPTreeNodeItem item) {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public boolean deleteItem(int index) {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public void printItems() {
        //To change body of implemented methods use File | Settings | File Templates.
    }
}
