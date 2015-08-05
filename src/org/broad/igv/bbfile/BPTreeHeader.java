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
import htsjdk.tribble.util.LittleEndianInputStream;

import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: martind
 * Date: Jan 8, 2010
 * Time: 3:50:08 PM
 * To change this template use File | Settings | File Templates.
 */
/*
 *  Container class for BBFile B+ Tree header.
 *  B+ Tree Header can be constructed by reading values in from a BBFile
  *  (  Table E), or by assigning the values in a constructor.
 *
 * */
public class  BPTreeHeader {

    private static Logger log = Logger.getLogger(BPTreeHeader.class);

    static public final int BPTREE_HEADER_SIZE = 32;

    static public final int BPTREE_MAGIC_LTH = 0x78CA8C91;
    static public final int BPTREE_MAGIC_HTL = 0x918CCA78;

    private long headerOffset;     // BBFile file offset for mChromosome tree
    private boolean headerOK;      // B+ Tree header OK?
    
    // Chromosome B+ Tree Header - Table E
    private int magic;        // magic number identifies it as B+ header
    private int blockSize;    // number of children per block
    private int keySize;      // min # of charcter bytes for mChromosome name
    private int valSize;      // size of (bigWig) values - currently 8
    private long itemCount;   // number of chromosomes/contigs in B+ tree
    private long reserved;    // Currently 0

   /*
   *    Constructor for reading in a B+ tree header a from a file input stream.
   *
   *    Parameters:
   *        fis - file input handle
   *        fileOffset - file offset to the B+ tree header
   *        isLowToHigh - indicates byte order is low to high, else is high to low
   * */
    public BPTreeHeader(SeekableStream fis, long fileOffset, boolean isLowToHigh) {

        long itemsCount;

       // save the seekable file handle  and B+ Tree file offset
       headerOffset = fileOffset;

       // Note: a bad B+ Tree header will result in false returned
       headerOK =  readHeader(fis, headerOffset, isLowToHigh);
    }



    public static int getHeaderSize() {
        return BPTREE_HEADER_SIZE;
    }

    public long getHeaderOffset() {
        return headerOffset;
    }

     public boolean isHeaderOK() {
        return headerOK;
    }

    public int getMagic() {
        return magic;
    }

    public int getBlockSize() {
        return blockSize;
    }

    public int getKeySize() {
        return keySize;
    }

    public int getValSize() {
        return valSize;
    }

    public long getItemCount() {
        return itemCount;
    }

    public long getReserved() {
        return reserved;
    }

    // prints out the B+ Tree Header
     public void print() {

        // Chromosome B+ Tree  Header - BBFile Table E
        if(headerOK)
            log.debug("B+ Tree Header was read from file location " + headerOffset);
        log.debug(" Magic ID =" + magic);
        log.debug(" Block size = " + blockSize);
        log.debug(" Key size = " + keySize);
        log.debug(" Indexed value size = " + valSize);
        log.debug(" Item Count = " + itemCount);
        log.debug(" Reserved = " + reserved);
    }
    
   /*
   * Reads in the B+ Tree Header.
   * Returns status of B+ tree header read; true if read, false if not.
   * */
    private boolean readHeader(SeekableStream fis, long fileOffset, boolean isLowToHigh) {

        LittleEndianInputStream lbdis;
        DataInputStream bdis;

         byte[] buffer = new byte[BPTREE_HEADER_SIZE];
         int bytesRead;
    
        try {
            // Read B+ tree header into a buffer
            fis.seek(fileOffset);
            fis.readFully(buffer);
        
            // decode header
            if(isLowToHigh){
                lbdis = new LittleEndianInputStream(new ByteArrayInputStream(buffer));

                // check for a valid B+ Tree Header
                magic = lbdis.readInt();

                if(magic != BPTREE_MAGIC_LTH)
                    return false;

                // Get mChromosome B+ header information
                blockSize = lbdis.readInt();
                keySize = lbdis.readInt();
                valSize = lbdis.readInt();
                itemCount = lbdis.readLong();
                reserved = lbdis.readLong();
            }
            else {
                bdis = new DataInputStream(new ByteArrayInputStream(buffer));

                // check for a valid B+ Tree Header
                magic = bdis.readInt();

                if(magic != BPTREE_MAGIC_HTL)
                    return false;

                // Get mChromosome B+ header information
                blockSize = bdis.readInt();
                keySize = bdis.readInt();
                valSize = bdis.readInt();
                itemCount = bdis.readLong();
                reserved = bdis.readLong();

            }

        }catch(IOException ex) {
           log.error("Error reading B+ tree header " + ex);
           throw new RuntimeException("Error reading B+ tree header \n", ex);
            }

        // success
         return true;
    }

}

