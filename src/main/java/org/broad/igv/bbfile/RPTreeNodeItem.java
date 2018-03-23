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
 * Time: 4:29:53 PM
 * To change this template use File | Settings | File Templates.
 */
/*
*   RPTreeNodeItem interface for storage of R+ tree node item information.
*
*   Note: The bounding 1D rectangle defined by:
*       (mStartBase mChromosome, mStartBase base) to (mEndBase mChromosome, mEndBase base)
*    is used as a key for insertion and searches on the R+ tree.
*
* */

abstract class RPTreeNodeItem {

    private static Logger log = Logger.getLogger(RPTreeNodeItem.class);


    protected RPChromosomeRegion chromosomeBounds; // chromosome bounds for item


    public RPTreeNodeItem(RPChromosomeRegion chromosomeBounds) {
        this.chromosomeBounds = chromosomeBounds;
    }

    public RPChromosomeRegion getChromosomeBounds() {
        return chromosomeBounds;
    }

    // Note: compareRegions returns the following values:
    //   -2 indicates chromosome region is completely below node region
    //   -1 indicates that chromosome region intersect node region from below
    //  0 means that chromosome region is inclusive to node region
    //  1 indicates chromosome region intersects node region from above
    //  2 indicates that this region is completely above that region

    public int compareRegions(RPChromosomeRegion chromosomeRegion) {
        return chromosomeBounds.compareRegions(chromosomeRegion);
    }

    public void print() {

        log.debug("Child node item :\n");
        log.debug(" StartChromID = " + chromosomeBounds.getStartChromID() + "\n");
        log.debug(" StartBase = " + chromosomeBounds.getStartBase() + "\n");
        log.debug(" EndChromID = " + chromosomeBounds.getEndChromID() + "\n");
        log.debug(" EndBase = " + chromosomeBounds.getEndBase() + "\n");


    }
}
