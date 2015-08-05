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

import htsjdk.samtools.seekablestream.SeekableStream;
import org.apache.log4j.Logger;

import java.util.*;

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
    private RPChromosomeRegion selectionRegion;  // selection region for iterator
    private boolean contained; // if true, features must be fully contained by extraction region

    // File access variables for reading Bed data block
    private SeekableStream fis;  // file input stream handle
    private BPTree chromIDTree;    // B+ chromosome index tree
    private RPTree chromDataTree;  // R+ chromosome data location tree

    // chromosome region extraction items
    private HashMap<Integer, String> chromosomeMap;  // map of chromosome ID's and corresponding names
    List<BedFeature> features;
    int currentIdx = 0;

    /**
     * Constructor for a BigBed iterator over the specified chromosome region
     * <p/>
     * Parameters:
     * fis - file input stream handle
     * chromIDTree - B+ index tree returns chromomosme ID's for chromosome names
     * chromDataTree - R+ chromosome data locations tree
     * selectionRegion - chromosome region for selection of Bed feature extraction
     * consists of:
     * startChromID - ID of start chromosome
     * startBase - starting base position for features
     * endChromID - ID of end chromosome
     * endBase - starting base position for features
     * contained - specifies bed features must be contained by region, if true;
     * else return any intersecting region features
     */
    public BigBedIterator(SeekableStream fis, BPTree chromIDTree, RPTree chromDataTree,
                          RPChromosomeRegion selectionRegion, boolean contained) {

        // check for valid selection region
        if (selectionRegion == null)
            throw new RuntimeException("Error: BigBedIterator selection region is null\n");

        this.fis = fis;
        this.chromIDTree = chromIDTree;
        this.chromDataTree = chromDataTree;
        this.selectionRegion = selectionRegion;
        this.contained = contained;

        List<RPTreeLeafNodeItem> leafNodeItems = chromDataTree.getChromosomeDataHits(selectionRegion, contained);
        features = new ArrayList<BedFeature>(512 * leafNodeItems.size());
        for (RPTreeLeafNodeItem item : leafNodeItems) {
            features.addAll(readBedDataBlock(item));
        }
    }

    public BigBedIterator() {
        features = Collections.emptyList();
    }


    public boolean hasNext() {
        return currentIdx < features.size();
    }

    public BedFeature next() {
        BedFeature retvalue = features.get(currentIdx);
        currentIdx++;
        return retvalue;
    }


    public void remove() {
        throw new UnsupportedOperationException("Remove iterator item is not supported yet.");
    }

    // ************ BigBedIterator specific methods *******************


    /*
    *   Method returns if bed items must be completely contained in
    *   the selection region.
    *
    *   Returns:
    *       Boolean indicates items must be contained in selection region if true,
    *       else may intersect the selection region if false
    * */
    public boolean isContained() {
        return contained;
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
    private List<BedFeature> readBedDataBlock(RPTreeLeafNodeItem leafHitItem) {

        // get the chromosome names associated with the hit region ID's
        int startChromID = leafHitItem.getChromosomeBounds().getStartChromID();
        int endChromID = leafHitItem.getChromosomeBounds().getEndChromID();
        chromosomeMap = chromIDTree.getChromosomeIDMap(startChromID, endChromID);

        boolean isLowToHigh = chromDataTree.isIsLowToHigh();
        int uncompressBufSize = chromDataTree.getUncompressBuffSize();

        // decompress leaf item data block for feature extraction
        BigBedDataBlock bedDataBlock = new BigBedDataBlock(fis, leafHitItem, chromosomeMap, isLowToHigh,
                uncompressBufSize);

        // get data block Bed feature list and set next index to first item
        return bedDataBlock.getBedData(selectionRegion, contained);

    }

}
