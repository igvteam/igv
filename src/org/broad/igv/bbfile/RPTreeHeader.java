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
import org.broad.tribble.util.LittleEndianInputStream;

import java.io.DataInputStream;
import java.io.IOException;
import java.io.ByteArrayInputStream;

/**
 * Created by IntelliJ IDEA.
 * User: martind
 * Date: Jan 14, 2010
 * Time: 3:59:55 PM
 * To change this template use File | Settings | File Templates.
 */
public class RPTreeHeader {

    private static Logger log = Logger.getLogger(RPTreeHeader.class);

    public final int RPTREE_HEADER_SIZE = 48;

    public final int RPTREE_MAGIC_LTH = 0x2468ACE0;
    public final int RPTREE_MAGIC_HTL = 0xE0AC6824;

    // defines the R+ Tree access

    private long rpTreeOffset;         // BBFile file offset for mChromosome region tree
    private boolean headerOK;          // R+ Tree header read OK

    // R+ Tree header - Table K
    private int magic;             // magic number identifies it as B+ header
    private int blockSize;         // number of children per block
    private long itemCount;        // number of chromosomes/contigs in B+ tree
    private int startChromID;      // ID of the first mChromosome in item
    private int startBase;         // Position of first base in item
    private int endChromID;        // ID of the first mChromosome in item
    private int endBase;           // Position of first base in item
    private long endFileOffset;    // file position marking mEndBase of data
    private int itemsPerSlot;      // number of items per leaf
    private long reserved;         // Currently 0

    // constructor   - reads from file input stream
    /*
    *   Constructor
    *
    *   Parameters:
    *       fis - file input stream handle
    *       fileOffset - file offset to the RP tree header
    *       isLowToHigh - if true, indicates low to high byte order, else high to low
    * */
    public RPTreeHeader(SeekableStream fis, long fileOffset, boolean isLowToHigh) {

        long itemsCount;

       rpTreeOffset =  fileOffset;

       // Note: a bad R+ Tree header will result in false returned
       headerOK =  readHeader(fis, rpTreeOffset, isLowToHigh);

    }

    public boolean isHeaderOK() {
        return headerOK;
    }

    public int getHeaderSize() {
        return RPTREE_HEADER_SIZE;
    }

    public long getTreeOffset() {
        return rpTreeOffset;
    }

    public int getMagic() {
        return magic;
    }

    public int getBlockSize() {
        return blockSize;
    }

    public long getItemCount() {
        return itemCount;
    }

    public int getStartChromID() {
        return startChromID;
    }

    public int getStartBase() {
        return startBase;
    }

    public int getEndChromID() {
        return endChromID;
    }

    public int getEndBase() {
        return endBase;
    }

    public long getMEndFileOffset() {
        return endFileOffset;
    }

    public int getItemsPerSlot() {
        return itemsPerSlot;
    }

    public long getReserved() {
        return reserved;
    }

// prints out the B+ Tree Header
public void print() {

   // note if read successfully
   if(headerOK){
       log.debug("R+ tree header has " + RPTREE_HEADER_SIZE + " bytes.");
       log.debug("R+ tree header magic = " + magic);
   }
   else {
       log.debug("R+ Tree header is unrecognized type, header magic = " + magic);
       return;
   }

   // Table E - Chromosome B+ Tree  Header
   log.debug("R+ Tree file offset = " + rpTreeOffset);
   log.debug("magic = " + magic);
   log.debug("Block size = " + blockSize);
   log.debug("ItemCount = " + itemCount);
   log.debug("StartChromID = " + startChromID);
   log.debug("StartBase = " + startBase);
   log.debug("EndChromID = " + endChromID);
   log.debug("EndBase = " + endBase);
   log.debug("EndFileOffset = " + endFileOffset);
   log.debug("ItemsPerSlot = " + itemsPerSlot);
   log.debug("Reserved = " + reserved);
   }

  /*
  * Reads in the R+ Tree Header.
  *
  * Returns status of for tree header read; true if read, false if not.
  * */
   private boolean readHeader(SeekableStream fis, long fileOffset, boolean isLowToHigh){

   LittleEndianInputStream lbdis;
   DataInputStream bdis;
      
    byte[] buffer = new byte[RPTREE_HEADER_SIZE];

   try {
       // Read R+ tree header into a buffer
       fis.seek(fileOffset);
       fis.readFully(buffer);

       // decode header
       if(isLowToHigh){
           lbdis = new LittleEndianInputStream(new ByteArrayInputStream(buffer));
           magic = lbdis.readInt();

           // check for a valid B+ Tree Header
           if(magic != RPTREE_MAGIC_LTH)
               return false;

           // Get mChromosome B+ header information
           blockSize = lbdis.readInt();
           itemCount = lbdis.readLong();
           startChromID = lbdis.readInt();
           startBase = lbdis.readInt();
           endChromID = lbdis.readInt();
           endBase = lbdis.readInt();
           endFileOffset = lbdis.readLong();
           itemsPerSlot = lbdis.readInt();
           reserved = lbdis.readInt();
       }
       else {
           bdis = new DataInputStream(new ByteArrayInputStream(buffer));

           // check for a valid B+ Tree Header
           magic = bdis.readInt();

           if(magic != RPTREE_MAGIC_HTL)
               return false;

           // Get mChromosome B+ header information
           blockSize = bdis.readInt();
           itemCount = bdis.readLong();
           startChromID = bdis.readInt();
           startBase = bdis.readInt();
           endChromID = bdis.readInt();
           endBase = bdis.readInt();
           endFileOffset = bdis.readLong();
           itemsPerSlot = bdis.readInt();
           reserved = bdis.readInt();
       }

   }catch(IOException ex) {
           log.error("Error reading R+ tree header " + ex);
           throw new RuntimeException("Error reading R+ tree header ", ex);
       }

   // success
    return true;
   }


}
