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
import org.broad.igv.util.CompressionUtils;
import org.broad.tribble.util.SeekableStream;

import java.io.IOException;
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
    //private SeekableStream fis;  // file input stream handle
    private BPTree chromIDTree;    // B+ chromosome index tree
    private RPTree chromDataTree;  // R+ chromosome data location tree

    // chromosome region extraction items
    private HashMap<Integer, String> chromosomeMap;  // map of chromosome ID's and corresponding names
    List<BedFeature> features;
    int currentIdx = 0;

    CompressionUtils compressionUtils;

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
    public BigBedIterator(String path, BPTree chromIDTree, RPTree chromDataTree,
                          RPChromosomeRegion selectionRegion, boolean contained, CompressionUtils compressionUtils){

        // check for valid selection region
        if (selectionRegion == null)
            throw new RuntimeException("Error: BigBedIterator selection region is null\n");

        this.compressionUtils = compressionUtils;

        this.chromIDTree = chromIDTree;
        this.chromDataTree = chromDataTree;
        this.selectionRegion = selectionRegion;
        this.contained = contained;

        List<RPTreeLeafNodeItem> leafNodeItems = chromDataTree.getChromosomeDataHits(selectionRegion, contained);
        features = new ArrayList<BedFeature>(512 * leafNodeItems.size());

        SeekableStream fis = null;
        try {
            fis = BBFileReader.getStream(path);
            for (RPTreeLeafNodeItem item : leafNodeItems) {
                features.addAll(readBedDataBlock(item, fis));
            }
        } catch (IOException e) {
            log.error(e);
            throw new RuntimeException(e);
        }finally {
            if(fis != null) try{
                fis.close();
            }catch (IOException e){
                log.error(e);
            }
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
    private List<BedFeature> readBedDataBlock(RPTreeLeafNodeItem leafHitItem, SeekableStream fis) {

        // get the chromosome names associated with the hit region ID's
        int startChromID = leafHitItem.getChromosomeBounds().getStartChromID();
        int endChromID = leafHitItem.getChromosomeBounds().getEndChromID();
        chromosomeMap = chromIDTree.getChromosomeIDMap(startChromID, endChromID);

        boolean isLowToHigh = chromDataTree.isIsLowToHigh();
        int uncompressBufSize = chromDataTree.getUncompressBuffSize();

        // decompress leaf item data block for feature extraction
        BigBedDataBlock bedDataBlock = new BigBedDataBlock(fis, leafHitItem, chromosomeMap, isLowToHigh,
                uncompressBufSize, compressionUtils);

        // get data block Bed feature list and set next index to first item
        return bedDataBlock.getBedData(selectionRegion, contained);

    }

}
