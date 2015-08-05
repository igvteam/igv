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

import java.util.ArrayList;
import java.io.IOException;
import java.util.HashMap;

/**
 * Created by IntelliJ IDEA.
 * User: martind
 * Date: Apr 16, 2010
 * Time: 10:32:51 AM
 * To change this template use File | Settings | File Templates.
 */
public class BigWigDataBlock {

    private static Logger log = Logger.getLogger(BigWigDataBlock.class);

    // BigWig data types sizes
    final int FIXED_STEP_ITEM_SIZE = 4;
    final int VAR_STEP_ITEM_SIZE = 8;
    final int BED_GRAPH_ITEM_SIZE = 12;

    // Bed data block access variables   - for reading in bed records from a file
    private long fileOffset;       // Wig data block file offset
    private long leafDataSize;     // byte size for data block specified in the R+ leaf
    private boolean isLowToHigh;   // if true, data is low to high byte order; else high to low

    // defines the bigWig data source
    private HashMap<Integer, String> chromosomeMap;  // map of chromosome ID's and corresponding names
    private RPTreeLeafNodeItem leafHitItem;   // R+ leaf item containing data block location

    // uncompressed byte stream buffer and readers
    private byte[] wigBuffer;      // buffer containing leaf block data uncompressed
    private int remDataSize;       // number of uncompressed data bytes not extracted

    // Wig data extraction members
    private ArrayList<WigItem> wigItemList;  // array of Wig section items

    /*
    *   Constructor for Wig data block reader.
    *
    *   Parameters:
    *       fis - file input stream handle
    *       leafHitItem - R+ tree leaf hit item containing data block file location and hit status
    *       chromIDTree - B+ chromosome index tree returns chromosome ID's for names
    *       isLowToHigh - byte order is low to high if true; else high to low
    *       uncompressBufSize - byte size for decompression buffer; else 0 for uncompressed
    *
    * */
    public BigWigDataBlock(SeekableStream fis, RPTreeLeafNodeItem leafHitItem,
                           HashMap<Integer, String> chromosomeMap, boolean isLowToHigh, int uncompressBufSize){
        this.leafHitItem = leafHitItem;
        this.chromosomeMap = chromosomeMap;
        this.isLowToHigh = isLowToHigh;

        fileOffset = this.leafHitItem.getDataOffset();
        leafDataSize = this.leafHitItem.geDataSize();
        byte[] buffer = new byte[(int) leafDataSize];

        // read Wig data block into a buffer
        try {
            fis.seek(fileOffset);
            fis.readFully(buffer);

            // decompress if necessary - the buffer size is 0 for uncompressed data
            // Note:  BBFile Table C specifies a decompression buffer size
            if(uncompressBufSize > 0)
                wigBuffer = (new CompressionUtils()).decompress(buffer, uncompressBufSize);
            else
                wigBuffer = buffer;    // use uncompressed read buffer directly
        }catch(IOException ex) {
            log.error("Error reading Wig section for leaf item ", ex);
            String error = String.format("Error reading Wig section for leaf item %d\n");
            throw new RuntimeException(error, ex);
        }

        // initialize unread data size
        remDataSize = wigBuffer.length;

        // use getWigData to extract data block items
    }

    /*
    *   Method reads all Wig data sections within the decompressed block buffer
    *   and returns those items in the chromosome selection region.
    *
    *   Parameters:
    *       selectionRegion - chromosome region for selecting Wig values
   *       contained - indicates selected data must be contained in selection region
    *           if true, else may intersect selection region
    *
    *   Returns:
    *      Wig sections in selected from the data block; else null for none selected.
    *
    * */
    public ArrayList<WigItem> getWigData(RPChromosomeRegion selectionRegion,
                                         boolean contained){

        wigItemList = new ArrayList<WigItem>();

        for(int index = 0; remDataSize > 0; ++index) {

            // extract items in the Wig data section
            // Note: A RuntimeException is thrown if wig section is not read properly
            BigWigSection wigSection = new BigWigSection(wigBuffer, chromosomeMap, isLowToHigh, leafHitItem);

            // get wig section items and section bytes read
            int sectionBytes = wigSection.getSectionData(selectionRegion, contained, wigItemList);

            // adjust remaining data block size
            remDataSize -= sectionBytes;
        }

        return wigItemList;
    }

    public void print() {

        log.debug("Wig section data referenced by leaf item ");

        for(int index = 0; index <= wigItemList.size(); ++index) {
            // BigWig sections print themselves
            wigItemList.get(index).print();
        }
    }


}
