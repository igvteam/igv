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
import org.broad.tribble.util.LittleEndianInputStream;
import org.broad.tribble.util.SeekableStream;

import java.util.ArrayList;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.DataInputStream;
import java.util.HashMap;


/**
 * Created by IntelliJ IDEA.
 * User: martind
 * Date: Jan 26, 2010
 * Time: 12:18:32 PM
 * To change this template use File | Settings | File Templates.
 */
/*
*   Container class for reading and storing a block of bed data items.
* */
public class BigBedDataBlock {

 private static Logger log = Logger.getLogger(BigBedDataBlock.class);

    // Bed data block access variables   - for reading in bed records from a file
    private SeekableStream fis;  // file input stream handle
    private long fileOffset;       // Bed data block file offset
    private long dataBlockSize;     // byte size for data block specified in the R+ leaf
    private boolean isLowToHigh;   // if true, data is low to high byte order; else high to low

    // defines the bigBed/bigWig source chromosomes
    private HashMap<Integer, String> chromosomeMap;  // map of chromosome ID's and corresponding names
    private RPTreeLeafNodeItem leafHitItem;   // R+ tree leaf item containing data block location

    // Provides uncompressed byte stream data reader
    private byte[] bedBuffer;  // buffer containing leaf block data uncompressed
    private int remDataSize;   // number of unread data bytes
    private long dataSizeRead;     // number of bytes read from the decompressed mWigBuffer

    // byte stream readers
    private LittleEndianInputStream lbdis;    // low to high byte stream reader
    private DataInputStream dis;       // high to low byte stream reader

    // Bed data extraction members
    private ArrayList<BedFeature> bedFeatureList; // array of BigBed data
    private int nItemsSelected;    // number of Bed features selected from this section

    /*
    *   Constructor for Bed data block reader.
    *
    *   Parameters:
    *       fis - file input stream handle
    *       leafItem - R+ tree leaf item containing chromosome region and file data location
    *       chromIDTree - B+ chromosome index tree returns chromosome ID's for names
    *       isLowToHigh - byte order is low to high if true; else high to low
    *       uncompressBufSize - byte size for decompression buffer; else 0 for uncompressed
    * */
    public BigBedDataBlock(SeekableStream fis, RPTreeLeafNodeItem leafHitItem,
            HashMap<Integer, String> chromosomeMap, boolean isLowToHigh, int uncompressBufSize){
        this.fis = fis;
        this.leafHitItem = leafHitItem;
        this.chromosomeMap = chromosomeMap;
        this.isLowToHigh = isLowToHigh;

        dataBlockSize = this.leafHitItem.geDataSize();
        byte[] buffer = new byte[(int) dataBlockSize];

        fileOffset = this.leafHitItem.getDataOffset();

        // read Bed data block into a buffer
        try {
            fis.seek(fileOffset);
            fis.readFully(buffer);

            // decompress if necessary - the buffer size is 0 for uncompressed data
            // Note:  BBFile Table C specifies a decompression buffer size
            if(uncompressBufSize > 0)
              bedBuffer = BBCompressionUtils.decompress(buffer, uncompressBufSize);
            else
              bedBuffer = buffer;    // use uncompressed read buffer directly

       }catch(IOException ex) {
            String error = String.format("Error reading Bed data for leaf item %d \n");
            log.error(error, ex);
            throw new RuntimeException(error, ex);
       }

        // wrap the bed buffer as an input stream
        if(this.isLowToHigh)
            lbdis = new LittleEndianInputStream(new ByteArrayInputStream(bedBuffer));
        else
            dis = new DataInputStream(new ByteArrayInputStream(bedBuffer));

        // initialize unread data size
        remDataSize = bedBuffer.length;

        // use methods getBedData or getNextFeature to extract block data
    }

    /*
    *   Method returns all Bed features within the decompressed block buffer
    *
    *   Parameters:
    *       selectionRegion - chromosome region for selecting Bed features
    *       contained - indicates selected data must be contained in selection region
    *           if true, else may intersect selection region
    *
    *   Returns:
    *      Bed feature items in the data block
    *
    *   Note: Remaining bytes to data block are used to determine end of reading
    *   since a zoom record count for the data block is not known.
    * */
    public ArrayList<BedFeature> getBedData(RPChromosomeRegion selectionRegion,
                                                boolean contained) {
        int itemNumber = 0;
        int chromID, chromStart, chromEnd;
        String restOfFields;
        int itemHitValue;

        // chromID + chromStart + chromEnd + rest 0 byte
        // 0 byte for "restOfFields" is always present for bed data
        int minItemSize = 3 * 4 + 1;

        // allocate the bed feature array list
        bedFeatureList = new ArrayList<BedFeature>();

        // check if all leaf items are selection hits
        RPChromosomeRegion itemRegion = new RPChromosomeRegion( leafHitItem.getChromosomeBounds());
        int leafHitValue = itemRegion.compareRegions(selectionRegion);
        
        try {
            for(int index = 0; remDataSize > 0; ++index) {
                itemNumber = index + 1;

                // read in BigBed item fields - BBFile Table I
                if(isLowToHigh){
                    chromID = lbdis.readInt();
                    chromStart= lbdis.readInt();
                    chromEnd = lbdis.readInt();
                    restOfFields = lbdis.readString();
                }
                else{
                    chromID = dis.readInt();
                    chromStart= dis.readInt();
                    chromEnd = dis.readInt();
                    restOfFields = dis.readUTF();
                }

                if(leafHitValue == 0) {     // contained leaf region items always added
                    String chromosome = chromosomeMap.get(chromID);
                    BedFeature bbItem = new BedFeature(itemNumber, chromosome,
                         chromStart, chromEnd, restOfFields);
                    bedFeatureList.add(bbItem);
                }
                else {                      // test for hit
                    itemRegion = new RPChromosomeRegion(chromID, chromStart, chromID, chromEnd);
                    itemHitValue = itemRegion.compareRegions(selectionRegion);

                    // abs(itemHitValue) == 1 for intersection; itemHitValue == 0 for contained
                    if(!contained && Math.abs(itemHitValue) < 2 ||
                            itemHitValue == 0) {
                        // add bed feature to item selection list
                        String chromosome = chromosomeMap.get(chromID);
                        BedFeature bbItem = new BedFeature(itemNumber, chromosome,
                             chromStart, chromEnd, restOfFields);
                        bedFeatureList.add(bbItem);
                    }
                }

                // compute data block remainder from size of item read
                // todo: check that restOfFields.length() does not also include 0 byte terminator
                remDataSize -= minItemSize + restOfFields.length();
            }

        }catch(IOException ex) {
            log.error("Read error for Bed data item " + itemNumber);

            // accept this as an end of block condition unless no items were read
            if(itemNumber == 1)
                throw new RuntimeException("Read error for Bed data item " + itemNumber);
        }

        return bedFeatureList;
    }

    public void print() {

        log.debug("BigBed data for " + bedFeatureList.size() + " items");

        for(int index = 0; index <= bedFeatureList.size(); ++index) {
            // BigBed data items print themselves
            bedFeatureList.get(index).print();
         }
    }
}
