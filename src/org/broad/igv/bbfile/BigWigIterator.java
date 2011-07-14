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
import org.apache.log4j.Logger;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: martind
 * Date: Apr 13, 2010
 * Time: 12:34:16 PM
 * To change this template use File | Settings | File Templates.
 */
public class BigWigIterator implements Iterator<WigItem> {

    private static Logger log = Logger.getLogger(BigWigIterator.class);

    boolean empty = false;

    //specification of chromosome selection region
    private RPChromosomeRegion selectionRegion;  // selection region for iterator
    private boolean isContained;     // if true, features must be fully contained by selection region

    // File access variables for reading Bed data block
    private SeekableStream fis;  // file input stream handle
    private BPTree chromIDTree;    // B+ chromosome index tree
    private RPTree chromDataTree;  // R+ chromosome data location tree

    // chromosome region extraction items
    private ArrayList<RPTreeLeafNodeItem> leafHitList; // array of leaf hits for selection region items
    private HashMap<Integer, String> chromosomeMap;  // map of chromosome ID's and corresponding names
    private int leafItemIndex;   // index of current leaf item being processed from leaf hit list
    RPTreeLeafNodeItem leafHitItem;   // leaf item being processed by next
    private RPChromosomeRegion hitRegion;  // hit selection region for iterator

    // current data block processing members
    private BigWigDataBlock wigDataBlock;  // Wig data block with Wig records decompressed
    private boolean dataBlockRead;  // indicates successful read of data block
    private ArrayList<WigItem> wigItemList; // array of selected Wig values
    private int wigItemIndex;      // index of next Wig data item from the list

    /**
     * Constructor for a BigWig iterator over the specified chromosome region
     * <p/>
     * Parameters:
     * fis - file input stream handle
     * chromIDTree - B+ chromosome index tree provides chromosome ID's for chromosome names
     * chromDataTree - R+ chromosome data locations tree
     * selectionRegion - chromosome region for selection of Wig feature extraction
     * consists of:
     * startChromID - ID of start chromosome
     * startBase - starting base position for values
     * endChromID - ID of end chromosome
     * endBase - ending base position for values
     * contained - specifies wig values must be contained by region, if true;
     * else return any intersecting region values
     */
    public BigWigIterator(SeekableStream fis, BPTree chromIDTree, RPTree chromDataTree,
                          RPChromosomeRegion selectionRegion, boolean contained) {

        // check for valid selection region
        if (selectionRegion == null)
            throw new RuntimeException("Error: BigWigIterator selection region is null\n");

        this.fis = fis;
        this.chromIDTree = chromIDTree;
        this.chromDataTree = chromDataTree;
        this.selectionRegion = new RPChromosomeRegion(selectionRegion);
        isContained = contained;

        // set up hit list and first data block read
        int hitCount = getHitRegion(selectionRegion, contained);
        if (hitCount == 0) {
            empty = true;
        }

        // Ready for next() data extraction
    }

    /**
     * Constructor for an "empty" iterator
     */
    public  BigWigIterator () {
        empty = true;

    }

    /*
   *  Method returns status on a "next item" being available.
   *
   *  Return:
   *      True if a "next item" exists; else false.
   *
   *  Note: If "next" method is called for a false condition,
   *      an NoSuchElementException will be thrown.
   * */

    public boolean hasNext() {

        if (empty) return false;

        // first check if current segment can be read for next Wig item
        if (wigItemIndex < wigItemList.size())
            return true;

            // need to fetch next data block
        else if (leafItemIndex < leafHitList.size())
            return true;

        else
            return false;
    }

    /**
     * Method returns the current Wig item and advances to the next Wig record.
     * <p/>
     * Returns:
     * Wig item for current BigWig data record.
     * <p/>
     * Note: If "next" method is called when a "next item" does not exist,
     * an NoSuchElementException will be thrown.
     */
    public WigItem next() {

        // return next Wig item in list
        if (wigItemIndex < wigItemList.size())
            return (wigItemList.get(wigItemIndex++));

            // attempt to get next leaf item data block
        else {
            int nHits = getHitRegion(selectionRegion, isContained);

            if (nHits > 0) {
                // Note: getDataBlock initializes bed feature index to 0
                return (wigItemList.get(wigItemIndex++)); // return 1st Data Block item
            } else {
                String result = String.format("Failed to find data for wig region (%d,%d,%d,%d)\n",
                        hitRegion.getStartChromID(), hitRegion.getStartBase(),
                        hitRegion.getEndChromID(), hitRegion.getEndBase());
                log.error(result);

                return null;
                //throw new NoSuchElementException(result);
            }
        }
    }

    public void remove() {
        throw new UnsupportedOperationException("Remove iterator item is not supported yet.");
    }

    // ************ BigBedIterator specific methods *******************

    /*
    *   Method returns the iterator selection region.
    * */

    public RPChromosomeRegion getSelectionRegion() {
        return selectionRegion;
    }

/*
    *   Method provides the iterator with a new selection region.
    *
    *   Parameters:
    *      selectionRegion - chromosome region for selection of Bed feature extraction
    *      consists of:
    *          startChromID - ID of start chromosome
    *          startBase - starting base position for features
    *          endChromID - ID of end chromosome
    *          endBase - starting base position for features
    *      contained - specifies bed features must be contained by region, if true;
    *          else return any intersecting region features
    *
    *   Returns:
    *       number of chromosome regions found in the selection region
    * */

    public int setSelectionRegion(RPChromosomeRegion selectionRegion,
                                  boolean contained) {
        this.selectionRegion = selectionRegion;
        isContained = contained;

        // set up hit list and first data block read
        leafHitList = null;    // Must nullify existing hit list first!
        int hitCount = getHitRegion(selectionRegion, contained);
        if (hitCount == 0)   // no hits - no point in fetching data
            throw new RuntimeException("No wig data found in the selection region");

        // Ready for next() data extraction

        return hitCount;
    }

    /*
    *   Method returns if bed items must be completely contained in
    *   the selection region.
    *
    *   Returns:
    *       Boolean indicates items must be contained in selection region if true,
    *       else may intersect the selection region if false
    * */

    public boolean isContained() {
        return isContained;
    }

    /*
    *   Method returns the BigBed file input stream handle.
    *
    *   Returns:
    *       File input stream handle
    * */

    public SeekableStream getBBFis() {
        return fis;
    }

    /*
    *   Method returns the B+ chromosome index tree used for identifying
    *   chromosome ID's used to specify R+ chromosome data locations.
    *
    *   Returns:
    *       B+ chromosome index tree
    * */

    public BPTree getChromosomeIDTree() {
        return chromIDTree;
    }

    /*
    *   Method returns the R+ chromosome data tree used for identifying
    *   chromosome data locations for the selection region.
    *
    *   Returns:
    *       R+ chromosome data locations tree
    * */

    public RPTree getChromosomeDataTree() {
        return chromDataTree;
    }

    /*
    *   Method finds the chromosome data hit items for the current hit selection region,
    *   and loads first hit data.
    *
    *   Parameters:
    *       hitRegion - selection region for extracting hit items
    *       contained - indicates hit items must contained in selection region if true;
    *       and if false, may intersect selection region
    *   Note: The selection region will be limited to accommodate  mMaxLeafHits; which terminates
    *       selection at the leaf node at which maxLeafHits is reached. Total number of selected
    *       items may exceed maxLeafHits, but only by the number of leaves in the cutoff leaf node.
    *
    *   Returns:
    *       number of R+ chromosome data hits
    * */

    private int getHitRegion(RPChromosomeRegion hitRegion, boolean contained) {

        int hitCount = 0;

        // check if new hit list is needed
        if (leafHitList == null) {
            hitCount = getHitList(hitRegion, contained);
            if (hitCount == 0)
                return 0;   // no hit data found
        } else {
            hitCount = leafHitList.size() - leafItemIndex;
            if (hitCount == 0)
                return 0;   // hit list exhausted
        }

        // Perform a block read for starting base of selection region - use first leaf hit
        dataBlockRead = getDataBlock(leafItemIndex++);

        // try next item - probably intersection issue
        // Note: recursive call until a block is valid or hit list exhuasted
        if (!dataBlockRead)
            hitCount = getHitRegion(hitRegion, contained);

        return hitCount;
    }

    /*
    *   Method finds the chromosome data tree hit items for the current hit selection region.
    *
    *   Parameters:
    *       hitRegion - selection region for extracting hit items
    *       contained - indicates hit items must contained in selection region if true;
    *       and if false, may intersect selection region
    *
    *   Note: The selection region will be limited to accommodate  mMaxLeafHits; which terminates
    *       selection at the leaf node at which maxLeafHits is reached. Total number of selected
    *       items may exceed maxLeafHits, but only by the number of leaves in the cutoff leaf node.
    *
    *   Returns:
    *       number of R+ chromosome data hits
    * */

    private int getHitList(RPChromosomeRegion hitRegion, boolean contained) {

        // hit list for hit region; subject to mMaxLeafHits limitation
        leafHitList = chromDataTree.getChromosomeDataHits(hitRegion, contained);

        // check if any leaf items were selected
        int nHits = leafHitList.size();
        if (nHits == 0)
            return 0;
        else
            leafItemIndex = 0;    // reset hit item index to start of list

        // find hit bounds from first and last hit items
        int startChromID = leafHitList.get(0).getChromosomeBounds().getStartChromID();
        int startBase = leafHitList.get(0).getChromosomeBounds().getStartBase();
        int endChromID = leafHitList.get(nHits - 1).getChromosomeBounds().getEndChromID();
        int endBase = leafHitList.get(nHits - 1).getChromosomeBounds().getEndBase();

        // save hit region; not currently used but useful for debug
        this.hitRegion = new RPChromosomeRegion(startChromID, startBase, endChromID, endBase);

        return nHits;
    }

    /*
    *   Method sets up a decompressed data block of big bed features for iteration.
    *
    *   Parameters:
    *       leafIteIndex - leaf item index in the hit list referencing the data block
    *
    *   Returns:
    *       Successful Bed feature data block set up: true or false.
    * */

    private boolean getDataBlock(int leafItemIndex) {

        // check for valid data block
        if (leafItemIndex >= leafHitList.size())
            return false;

        // Perform a block read for indexed leaf item
        leafHitItem = leafHitList.get(leafItemIndex);

        // get the chromosome names associated with the hit region ID's
        int startChromID = leafHitItem.getChromosomeBounds().getStartChromID();
        int endChromID = leafHitItem.getChromosomeBounds().getEndChromID();
        chromosomeMap = chromIDTree.getChromosomeIDMap(startChromID, endChromID);

        boolean isLowToHigh = chromDataTree.isIsLowToHigh();
        int uncompressBufSize = chromDataTree.getUncompressBuffSize();

        // decompress leaf item data block for feature extraction
        wigDataBlock = new BigWigDataBlock(fis, leafHitItem, chromosomeMap, isLowToHigh, uncompressBufSize);

        // get section Wig item list and set next index to first item
        wigItemList = wigDataBlock.getWigData(selectionRegion, isContained);
        wigItemIndex = 0;

        // data block items available for iterator
        if (wigItemList.size() > 0)
            return true;
        else
            return false;
    }

}
