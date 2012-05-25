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

import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by IntelliJ IDEA.
 * User: Owner
 * Date: May 3, 2010
 * Time: 12:32:17 AM
 * To change this template use File | Settings | File Templates.
 */
public class BigWigSection {

    private static Logger log = Logger.getLogger(BigWigSection.class);

    private boolean isLowToHigh;       // byte order is low to high if true; else high to low
    private LittleEndianInputStream lbdis;   // input stream reader for low to high byte ordered data
    private DataInputStream dis;      // input stream reader for high to low byte ordered data

    private RPTreeLeafNodeItem leafHitItem;    // leaf item defines chromosome region and file data location
    private int sectionDataSize;       // byte size of decompressed data for this section
    private HashMap<Integer, String> chromosomeMap; // map of chromosome ID's and corresponding names
    private BigWigSectionHeader wigSectionHeader;  // wig section header

    /*
    *   Constructor for a BigWig data section which includes the section header
    *   and Wig data items.
    *
    *   Parameters:
    *       sectionBuffer - buffer contains decompressed Wig section header + data
    *       sectionIndex - wig section index for leaf data block
    *       chromIDTree - B+ chromosome index tree returns chromosome names for ID's
    *       isLowToHigh - if true, data byte order is low to high ; else is high to low
    *       leafHitItem - contains leaf node information for testing against selection region
    *
    * */
    public BigWigSection(byte[] sectionBuffer, HashMap<Integer, String> chromosomeMap,
                         boolean isLowToHigh, RPTreeLeafNodeItem leafHitItem){

        this.chromosomeMap =  chromosomeMap;
        this.isLowToHigh = isLowToHigh;
        this.leafHitItem = leafHitItem;

        // wrap the Wig section buffer as an input stream and get the section header
        // Note: A RuntimeException is thrown if header is not read properly
        if(this.isLowToHigh) {
            lbdis = new LittleEndianInputStream(new ByteArrayInputStream(sectionBuffer));
            wigSectionHeader = new BigWigSectionHeader(lbdis);

        }
        else  {
            dis = new DataInputStream(new ByteArrayInputStream(sectionBuffer));
            wigSectionHeader = new BigWigSectionHeader(dis);
        }

        // check for valid Wig item type
        if(wigSectionHeader.getItemType() == BigWigSectionHeader.WigItemType.Unknown){
            throw new RuntimeException("Read error on wig section leaf index ");
        }

        // include header in data segment size accounting
        sectionDataSize = wigSectionHeader.SECTION_HEADER_SIZE;

        // use method getSectionData to extract section data
    }

    /*
    *   Method returns the if the Wig items defined in this section are valid.
    *
    *   Note: Use BigWigSectionHeader to obtain more information on
    *       Wig section data specifications.
    *
    *   Returns:
    *       Specifies if Wig section has a valid data item type.
    * */
    public boolean isValidSectionType(){
        return wigSectionHeader.IsValidType();
    }

    /*
    *   Method returns the Wig Section Header
    *
    *   Returns:
    *       Wig section header
    * */
    public int getItemCount() {
        return wigSectionHeader.getItemCount();
    }

    /*
    *   Method returns the Wig Section Header
    *
    *   Returns:
    *       Wig section header
    * */
    public BigWigSectionHeader getSectionHeader() {
        return wigSectionHeader;
    }

    /*
    *   Method returns the number bytes of decompressed data in this section.
    *
    *   Returns:
    *       Number of uncompressed bytes read for the Wig data section
    * */
    public int getSectionDataSize() {
        return sectionDataSize;
    }

    /*
    *   Method reads Wig data items within the decompressed block buffer for the selection region.
    *
    *   Parameters:
    *       selectionRegion - chromosome selection region for item extraction
    *       contained - indicates select region must be contained in value region
    *           if true, else may intersect selection region for extraction
    *
    *   Returns:
    *     Size in bytes for the wig data section.
    *     Items read in the wig segment data block are added to the wig item list .
    *
     *   Note: Unlike ZoomLevel and BigBed formats, the Wig Section data block header contains
    *   an item count used to determine the end of data read.
    * */
    public int getSectionData(RPChromosomeRegion selectionRegion, boolean contained,
                              ArrayList<WigItem> wigItemList) {

        // get the section's data item specifications
        // Note: A RuntimeException is thrown if wig section header is not read properly
        int chromID =  wigSectionHeader.getChromID();
        String chromosome = chromosomeMap.get(chromID);
        int itemCount = wigSectionHeader.getItemCount();
        int chromStart = wigSectionHeader.getChromosomeStart();
        int chromEnd = wigSectionHeader.getChromosomeEnd();
        int itemStep = wigSectionHeader.getItemStep();
        int itemSpan =  wigSectionHeader.getItemSpan();
        int itemIndex = 0;
        int startBase = 0;
        int endBase = 0;
        float value = 0.0f;

        // find Wig data type - BBFile Table J item type
        BigWigSectionHeader.WigItemType itemType = wigSectionHeader.getItemType();

        // check if all leaf items are selection hits
        RPChromosomeRegion itemRegion = new RPChromosomeRegion(chromID, chromStart,
                            chromID, chromEnd);
        int leafHitValue = itemRegion.compareRegions(selectionRegion);


        // extract Wig data records
        // Note: the buffer input stream is positioned past section header
        try {
            for(int index = 0; index < itemCount; ++index) {
                ++itemIndex;
                if(isLowToHigh){
                    if(itemType == BigWigSectionHeader.WigItemType.FixedStep){
                        startBase = chromStart;
                        endBase = startBase + itemSpan;
                        value = lbdis.readFloat();
                        chromStart = startBase + itemStep;
                        sectionDataSize += BigWigSectionHeader.FIXEDSTEP_ITEM_SIZE;
                    }
                    else if(itemType == BigWigSectionHeader.WigItemType.VarStep){

                        startBase = lbdis.readInt();
                        endBase = startBase + itemSpan;
                        value = lbdis.readFloat();
                        sectionDataSize += BigWigSectionHeader.VARSTEP_ITEM_SIZE;
                    }
                    else if(itemType == BigWigSectionHeader.WigItemType.BedGraph){
                        startBase = lbdis.readInt();
                        endBase = lbdis.readInt();
                        value = lbdis.readFloat();
                        sectionDataSize += BigWigSectionHeader.BEDGRAPH_ITEM_SIZE;
                    }
                }
                else {  // byte order is high to low
                    if(itemType == BigWigSectionHeader.WigItemType.FixedStep){
                        startBase = chromStart;
                        endBase = startBase + itemSpan;
                        value = dis.readFloat();
                        chromStart = startBase + itemStep;
                        sectionDataSize += BigWigSectionHeader.FIXEDSTEP_ITEM_SIZE;
                    }
                    else if(itemType == BigWigSectionHeader.WigItemType.VarStep){
                        startBase = dis.readInt();
                        endBase = startBase + itemSpan;
                        value = dis.readFloat();
                        sectionDataSize += BigWigSectionHeader.VARSTEP_ITEM_SIZE;
                    }
                    else if(itemType == BigWigSectionHeader.WigItemType.BedGraph){
                        startBase = dis.readInt();
                        endBase = dis.readInt();
                        value = dis.readFloat();
                        sectionDataSize += BigWigSectionHeader.BEDGRAPH_ITEM_SIZE;
                    }
                }

                // contained leaf region items are always added - otherwise test conditions
                if(leafHitValue == 0) {
                    WigItem bbItem = new WigItem(itemIndex, chromosome, startBase, endBase, value);
                    wigItemList.add(bbItem);
                }
                else {
                    itemRegion = new RPChromosomeRegion(chromID, startBase, chromID, endBase);
                    int itemHitValue = itemRegion.compareRegions(selectionRegion);

                    // hitValue < 2 needed for intersection; hitValue < 1 needed for contained = true
                    if(itemHitValue == 0 || !contained && Math.abs(itemHitValue) < 2) {
                        WigItem bbItem = new WigItem(itemIndex, chromosome, startBase, endBase, value);
                        wigItemList.add(bbItem);
                    }
                }

            }

        }catch(IOException ex) {
            log.error("Read error for Wig section item " + itemIndex);
            throw new RuntimeException("Read error for Wig section item " + itemIndex);
        }

        return sectionDataSize;
    }

    /*
    *   Method prints out the data items for this Wig section.
    * */
    public void print() {
        log.debug("Wig section for leaf item  has a data size = " + sectionDataSize);
        wigSectionHeader.print();
    }

}
