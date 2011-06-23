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
import org.broad.tribble.util.SeekableStream;

import java.util.HashMap;
import java.util.Iterator;
import java.util.ArrayList;
import java.util.NoSuchElementException;

/**
 * Created by IntelliJ IDEA.
 * User: martind
 * Date: Apr 5, 2010
 * Time: 3:10:26 PM
 * To change this template use File | Settings | File Templates.
 */
public class BigBedIterator implements Iterator<BedFeature> {

     private static Logger log = Logger.getLogger(BigBedIterator.class);

    //specification of chromosome selection region
    private RPChromosomeRegion mSelectionRegion;  // selection region for iterator
    private boolean mIsContained; // if true, features must be fully contained by extraction region
    private RPChromosomeRegion mHitRegion;  // hit selection region for iterator

    // File access variables for reading Bed data block
    private SeekableStream mBBFis;  // file input stream handle
    private BPTree mChromIDTree;    // B+ chromosome index tree
    private RPTree mChromDataTree;  // R+ chromosome data location tree

    // chromosome region extraction items
    private ArrayList<RPTreeLeafNodeItem> mLeafHitList; // array of leaf hits for selection region items
    private HashMap<Integer, String> mChromosomeMap;  // map of chromosome ID's and corresponding names
    private int mLeafItemIndex;  // index of current leaf item being processed from leaf hit list
    RPTreeLeafNodeItem mLeafHitItem;   // leaf item being processed by next

    // current data block processing members
    private BigBedDataBlock mBedDataBlock; // Bed data block with Bed records decompressed
    private boolean mDataBlockRead;  // flag indicates successful read of data block
    private ArrayList<BedFeature> mBedFeatureList; // array of selected  Bed features
    private int mBedFeatureIndex;       // index of next Bed feature from the list

    /**
     *  Constructor for a BigBed iterator over the specified chromosome region
     *
     * Parameters:
     *      fis - file input stream handle
     *      chromIDTree - B+ index tree returns chromomosme ID's for chromosome names
     *      chromDataTree - R+ chromosome data locations tree
     *      selectionRegion - chromosome region for selection of Bed feature extraction
     *      consists of:
     *          startChromID - ID of start chromosome
     *          startBase - starting base position for features
     *          endChromID - ID of end chromosome
     *          endBase - starting base position for features
     *      contained - specifies bed features must be contained by region, if true;
     *          else return any intersecting region features
     */
     public BigBedIterator(SeekableStream fis, BPTree chromIDTree, RPTree chromDataTree,
            RPChromosomeRegion selectionRegion, boolean contained) {

        // check for valid selection region
        if(selectionRegion == null)
            throw new RuntimeException("Error: BigBedIterator selection region is null\n");

        mBBFis = fis;
        mChromIDTree = chromIDTree;
        mChromDataTree = chromDataTree;
        mSelectionRegion = selectionRegion;
        mIsContained = contained;

        // set up hit list and first data block read
        int hitCount = getHitRegion(selectionRegion, contained);
        if(hitCount == 0)   // no hits - no point in fetching data
            throw new RuntimeException("No wig data found in the selection region");

        // Ready for next() data extraction

    }

     /*
     *  Method returns status on a "next item" being available.
     *
     *  Return:
     *      True if a "next item" exists; else false.
     *
     *  Note: If "next" method is called for a false condition,
     *      an UnsupportedOperationException will be thrown.
     * */
     public boolean hasNext() {

        // first check if current data block can be read for next
        if(mBedFeatureIndex < mBedFeatureList.size())
            return true;

        // need to fetch next data block
        else if(mLeafItemIndex < mLeafHitList.size())
            return true;
 
         else
            return false;
    }

    /**
     *  Method returns the current bed feature and advances to the next bed record.
     *
     *  Returns:
     *      Bed feature for current BigBed data record.
     *
     *  Note: If "next" method is called when a "next item" does not exist,
     *      an UnsupportedOperationException will be thrown.
    */
    public BedFeature next() {

        // Is there a need to fetch next data block?
        if(mBedFeatureIndex < mBedFeatureList.size())
            return(mBedFeatureList.get(mBedFeatureIndex++));

        // attempt to get next leaf item data block
        else {
            int nHits = getHitRegion(mSelectionRegion, mIsContained);

            if(nHits > 0){
                // Note: getDataBlock initializes bed feature index to 0
                return(mBedFeatureList.get(mBedFeatureIndex++)); // return 1st Data Block item
            }
            else{
                 String result = String.format("Failed to find data for bed region (%d,%d,%d,%d)\n",
                    mHitRegion.getStartChromID(),  mHitRegion.getStartBase(),
                        mHitRegion.getEndChromID(), mHitRegion.getEndBase());
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
        return mSelectionRegion;
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
        mSelectionRegion = selectionRegion;
        mIsContained = contained;

        // set up hit list and first data block read
        mLeafHitList = null;    // Must nullify existing hit list first!
        int hitCount = getHitRegion(selectionRegion, contained);
        if(hitCount == 0)   // no hits - no point in fetching data
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
        return mIsContained;
    }

    /*
    *   Method returns the BigBed file input stream handle.
    *
    *   Returns:
    *       File input stream handle
    * */
    public SeekableStream getBBFis() {
        return mBBFis;
    }

    /*
    *   Method returns the B+ chromosome index tree used for identifying
    *   chromosome ID's used to specify R+ chromosome data locations.
    *
    *   Returns:
    *       B+ chromosome index tree
    * */
    public BPTree getChromosomeIDTree() {
        return mChromIDTree;
    }

    /*
    *   Method returns the R+ chromosome data tree used for identifying
    *   chromosome data locations for the selection region.
    *
    *   Returns:
    *       R+ chromosome data locations tree
    * */
    public RPTree getChromosomeDataTree() {
        return mChromDataTree;
    }

    /*
    *   Method returns leaf items for the chromosome selection region.
    *
    *   Returns:
    *       List of leaf items with data locations for the selection region.
    * */
    public ArrayList<RPTreeLeafNodeItem> getLeafItems() {
        return mLeafHitList;
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
        if(mLeafHitList == null){
            hitCount = getHitList(hitRegion, contained);
            if(hitCount == 0)
                return 0;   // no hit data found
        }
         else {
            hitCount =  mLeafHitList.size() - mLeafItemIndex;
            if(hitCount == 0)
                return 0;   // hit list exhausted
        }

        // Perform a block read for starting base of selection region - use first leaf hit
        mDataBlockRead = getDataBlock(mLeafItemIndex++);

        // try next item - probably intersection issue
        // Note: recursive call until a block is valid or hit list exhuasted
        if(!mDataBlockRead)
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
         mLeafHitList = mChromDataTree.getChromosomeDataHits(hitRegion, contained);

        // check if any leaf items were selected
        int nHits = mLeafHitList.size();
        if(nHits == 0)
            return 0;   // no data hits found
        else
             mLeafItemIndex = 0;    // reset hit item index to start of list

        // find hit bounds
        int startChromID = mLeafHitList.get(0).getChromosomeBounds().getStartChromID();
        int startBase = mLeafHitList.get(0).getChromosomeBounds().getStartBase();
        int endChromID = mLeafHitList.get(nHits-1).getChromosomeBounds().getEndChromID();
        int endBase = mLeafHitList.get(nHits-1).getChromosomeBounds().getEndBase();

        // save hit region definition; not currently used but useful for debug
        mHitRegion = new  RPChromosomeRegion(startChromID, startBase, endChromID, endBase);

        return nHits;
    }

     /*
    *   Method sets up a decompressed data block of big bed features for iteration.
    *
    *   Parameters:
    *       leafItemIndex - leaf item index in the hit list referencing the data block
    *
    *   Returns:
    *       Successful Bed feature data block set up: true or false.
    * */
    private boolean getDataBlock(int leafItemIndex){

        // check for valid data block
        if(leafItemIndex >= mLeafHitList.size())
                return false;

        // Perform a block read for indexed leaf item
        mLeafHitItem = mLeafHitList.get(leafItemIndex);

        // get the chromosome names associated with the hit region ID's
        int startChromID = mLeafHitItem.getChromosomeBounds().getStartChromID();
        int endChromID = mLeafHitItem.getChromosomeBounds().getEndChromID();
        mChromosomeMap = mChromIDTree.getChromosomeIDMap(startChromID, endChromID);

        boolean isLowToHigh = mChromDataTree.isIsLowToHigh();
        int uncompressBufSize = mChromDataTree.getUncompressBuffSize();

        // decompress leaf item data block for feature extraction
        mBedDataBlock = new BigBedDataBlock(mBBFis, mLeafHitItem, mChromosomeMap, isLowToHigh,
                uncompressBufSize);

        // get data block Bed feature list and set next index to first item
        mBedFeatureList =  mBedDataBlock.getBedData(mSelectionRegion, mIsContained);
        mBedFeatureIndex = 0;

         // data block items available for iterator
        if(mBedFeatureList.size() > 0)
            return true;
        else
            return false;
    }

}
