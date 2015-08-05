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
import htsjdk.tribble.util.LittleEndianInputStream;

import java.io.IOException;
import java.io.DataInputStream;

/**
 * Created by IntelliJ IDEA.
 * User: martind
 * Date: Jan 6, 2010
 * Time: 3:00:11 PM
 * To change this template use File | Settings | File Templates.
 */
/*
*   Container class for BigWig section header class for data items - BBFile Table J
*
*   Note: appropriate WIG data formats are accomodated
*   according to WIG type in Table J
* */
public class BigWigSectionHeader {

    public enum WigItemType {
        BedGraph,
        VarStep,
        FixedStep,
        Unknown     // bad value
    }

    private static Logger log = Logger.getLogger(BigWigSectionHeader.class);

    public static final int SECTION_HEADER_SIZE = 24;
    public static final int FIXEDSTEP_ITEM_SIZE = 4;
    public static final int VARSTEP_ITEM_SIZE = 8;
    public static final int BEDGRAPH_ITEM_SIZE = 12;

    private int chromID;       // Chromosome/contig Numerical ID from BBFile Chromosome B+ tree
    private int chromStart;    // starting base position
    private int chromEnd;      // ending base position
    private int itemStep;      // number of base spaces between fixed items
    private int itemSpan;      // number of bases in fixed step items
    private WigItemType itemType; // type of data items: 1 = bedGraph, 2 = varStep, 3 = fixedStep
    private byte reserved;     // reserved; currently = 0
    private short itemCount;   // number of data items in this chromosome section

    private boolean isValidType;    // indicates a if a valid Wig item type was read
    private String itemDescription; // string representation of item type.

    /*
    *   Constructor creates a Wig Section Header (Table J) from uncompressed buffer.
    *
    *   Parameters:
    *       mLbdis - // buffer stream containing section header arranged low to high bytes
    * */
    public BigWigSectionHeader(LittleEndianInputStream lbdis) {

        byte type;

        // get Wig Section Header
        try {
            chromID = lbdis.readInt();
            chromStart = lbdis.readInt();
            chromEnd = lbdis.readInt();
            itemStep = lbdis.readInt();
            itemSpan = lbdis.readInt();
            type = lbdis.readByte();
            reserved = lbdis.readByte();
            itemCount = lbdis.readShort();
        }catch(IOException ex) {
            log.error("Error reading wig section header ", ex);
            throw new RuntimeException("Error reading wig section header", ex);
        }

        // tag as valid
        isValidType = getItemType(type);
    }

    /*
    *   Constructor creates a Wig Section Header (Table J) from uncompressed buffer.
    *
    *   Parameters:
    *       mLbdis - // buffer stream containing section header arranged high to low bytes
    * */
    public BigWigSectionHeader(DataInputStream bdis) {

        byte type;

        // get Wig Section Header
        try {
            chromID = bdis.readInt();
            chromStart = bdis.readInt();
            chromEnd = bdis.readInt();
            itemStep = bdis.readInt();
            itemSpan = bdis.readInt();
            type = bdis.readByte();
            reserved = bdis.readByte();
            itemCount = bdis.readShort();
        }catch(IOException ex) {
            log.error("Error reading wig section header ", ex);
            throw new RuntimeException("Error reading wig section header", ex);
        }

        // tag as valid
        isValidType = getItemType(type);
    }

    /*
    *   Method returns the chromosome ID
    *
    *   Returns:
    *       Chromosome ID for the section's region
    * */
    public int getChromID() {
        return chromID;
    }

    /*
    *   Method returns the chromosome starting base
    *
    *   Returns:
    *       Chromosome start base for the section's region
    * */
    public int getChromosomeStart() {
        return chromStart;
    }

    /*
    *   Method returns the chromosome ending base
    *
    *   Returns:
    *       Chromosome end base for the section's region
    * */
    public int getChromosomeEnd() {
        return chromEnd;
    }

    /*
    *   Method returns the base pairs step between items.
    *
    *   Returns:
    *       Chromosome base step between fixed step sections
    * */
    public int getItemStep() {
        return itemStep;
    }

    /*
    *   Method returns the base pairs span in items.
    *
    *   Returns:
    *       Chromosome base span for fixed and variable step sections
    * */
    public int getItemSpan() {
        return itemSpan;
    }

    /*
    *   Method returns the item type for the section's Wig data.
    *
    *   Returns:
    *       Section item type for Wig data
    * */
    public WigItemType getItemType() {
        return itemType;
    }

    /*
    *   Method returns if the section's data item type is valid.
    *
    *   Returns:
    *       Specifies if section's data iytem type is valid
    * */
    public boolean IsValidType() {
        return isValidType;
    }

    /*
    *   Method returns the number of section items.
    *
    *   Returns:
    *       Number of items defined for the section
    * */
    public short getItemCount() {
        return itemCount;
    }

    /*
    *   Method returns the reserved value for the section.
    *
    *   Returns:
    *       Reserved byte for the section (should always be 0)
    * */
    public byte getReserved() {
        return reserved;
    }

    public void print(){
        log.debug(" BigWig section header "
                + " for "+ itemDescription + " data");
        log.debug("Chromosome ID = " + chromID);
        log.debug("ChromStart = " + chromStart);
        log.debug("ChromEnd = " + chromEnd);
        log.debug("ItemStep = " + itemStep);
        log.debug("ItemSpan = " + itemSpan);
        log.debug("ItemType = " + itemType);
        log.debug("mReserved = " + reserved);
        log.debug("mItemCount = " + itemCount);
    }

    /*
    *   Method determines the Wig data type.
    *
    *   Parameters:
    *       byte type read from Wig section header
    *
    *   Returns:
    *       Indicates if type is a valid Wig item type
    * */
    private boolean getItemType(byte type){
        boolean isValid;

        if(type == 1){
            itemType = WigItemType.BedGraph;
            itemDescription = "Wig Bed Graph";
            isValid = true;
        }
        else if(type == 2){
            itemType = WigItemType.VarStep;
            itemDescription = "Wig Variable Step";
            isValid = true;
        }
        else if(type == 3){
            itemType = WigItemType.FixedStep;
            itemDescription = "Wig Fixed Step";
            isValid = true;
        }
        else {
            itemType = WigItemType.Unknown;
            itemDescription = "Wig Type Unknown";
            isValid = false;
        }

        return isValid;
    }

}
