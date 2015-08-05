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
import org.broad.igv.util.CompressionUtils;
import htsjdk.tribble.util.LittleEndianInputStream;

import java.io.ByteArrayOutputStream;
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
    private long fileOffset;       // Bed data block file offset
    private long dataBlockSize;     // byte size for data block specified in the R+ leaf
    private boolean isLowToHigh;   // if true, data is low to high byte order; else high to low

    // defines the bigBed/bigWig source chromosomes
    private HashMap<Integer, String> chromosomeMap;  // map of chromosome ID's and corresponding names
    private RPTreeLeafNodeItem leafHitItem;   // R+ tree leaf item containing data block location

    // Provides uncompressed byte stream data reader
    private byte[] bedBuffer;  // buffer containing leaf block data uncompressed
    private int remDataSize;   // number of unread data bytes

    // byte stream readers
    private LittleEndianInputStream lbdis;    // low to high byte stream reader
    private DataInputStream dis;       // high to low byte stream reader

    // Bed data extraction members
    private ArrayList<BedFeature> bedFeatureList; // array of BigBed data

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
                           HashMap<Integer, String> chromosomeMap, boolean isLowToHigh, int uncompressBufSize) {

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
            if (uncompressBufSize > 0)
                bedBuffer = (new CompressionUtils()).decompress(buffer, uncompressBufSize);
            else
                bedBuffer = buffer;    // use uncompressed read buffer directly

        } catch (IOException ex) {
            String error = String.format("Error reading Bed data for leaf item %d \n");
            log.error(error, ex);
            throw new RuntimeException(error, ex);
        }

        // wrap the bed buffer as an input stream
        if (this.isLowToHigh)
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

        // chromID + chromStart + chromEnd + rest 0 byte
        // 0 byte for "restOfFields" is always present for bed data
        int minItemSize = 3 * 4 + 1;

        // allocate the bed feature array list
        bedFeatureList = new ArrayList<BedFeature>();

        try {
            for (int index = 0; remDataSize >= minItemSize; ++index) {
                itemNumber = index + 1;

                // read in BigBed item fields - BBFile Table I
                if (isLowToHigh) {
                    chromID = lbdis.readInt();
                    chromStart = lbdis.readInt();
                    chromEnd = lbdis.readInt();
                    restOfFields = lbdis.readString();
                } else {
                    chromID = dis.readInt();
                    chromStart = dis.readInt();
                    chromEnd = dis.readInt();
                    restOfFields = readHighToLowString();
                }

                int leafHitValue = selectionRegion.compareRegions(chromID, chromStart, chromID, chromEnd);
                if (leafHitValue == -2 || (contained && leafHitValue == -1)) {
                    remDataSize -= (minItemSize + restOfFields.length());
                    break;
                }
                if (leafHitValue == 2 || (contained && leafHitValue == 1)) {
                    remDataSize -= (minItemSize + restOfFields.length());
                    continue;
                } else {
                    remDataSize -= (minItemSize + restOfFields.length());
                    String chromosome = chromosomeMap.get(chromID);
                    BedFeature bbItem = new BedFeature(itemNumber, chromosome, chromStart, chromEnd, restOfFields);
                    bedFeatureList.add(bbItem);
                }
            }

        } catch (IOException ex) {

            log.error("Read error for Bed data item " + itemNumber, ex);

            // accept this as an end of block condition unless no items were read
            if (itemNumber == 1)
                throw new RuntimeException("Read error for Bed data item " + itemNumber, ex);
        }

        return bedFeatureList;
    }

    private String readHighToLowString() throws IOException {
        String restOfFields;
        ByteArrayOutputStream bis = new ByteArrayOutputStream(100);
        byte b;
        while ((b = (byte) dis.read()) != 0) {
            bis.write(b);
        }
        restOfFields = new String(bis.toByteArray());
        return restOfFields;
    }

    public void print() {

        log.debug("BigBed data for " + bedFeatureList.size() + " items");

        for (int index = 0; index <= bedFeatureList.size(); ++index) {
            // BigBed data items print themselves
            bedFeatureList.get(index).print();
        }
    }
}
