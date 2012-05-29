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

/**
 * Created by IntelliJ IDEA.
 * User: martind
 * Date: Nov 20, 2009
 * Time: 3:49:14 PM
 * To change this template use File | Settings | File Templates.
 */

/**
 *  Container class defines the header information for BigBed and BigWig files
 */
package org.broad.igv.bbfile;


import org.apache.log4j.Logger;
import org.broad.tribble.util.SeekableStream;
import org.broad.tribble.util.*;
import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.IOException;

//import net.sf.samtools.util.SeekableStream;

/*
*   Container class for holding the BBFile header information, Table C .
**/
public class BBFileHeader {

    private static Logger log = Logger.getLogger(BBFileHeader.class);

    // defines bigBed/bigwig Header Format types
    static public final int BBFILE_HEADER_SIZE = 64;

    static public final int BIGWIG_MAGIC_LTH = 0x888FFC26; // BigWig Magic Low to High
    static public final int BIGWIG_MAGIC_HTL = 0x26FC8F66; // BigWig Magic High to Low

    static public final int BIGBED_MAGIC_LTH = 0x8789F2EB; // BigBed Magic Low to High
    static public final int BIGBED_MAGIC_HTL = 0xEBF28987; // BigBed Magic High to Low

    // defines the bigBed/bigWig source file access
    private String path;               // bigBed file/pathname
    private SeekableStream fis;      // BBFile I/O stream handle
    private long fileHeaderOffset;     // file offset for file header

    private boolean isHeaderOK;        // File header read correctly?
    private boolean isLowToHigh;       // flag indicates values represented low to high bytes
    private boolean isBigBed;          // flag indicates file is BigBed format
    private boolean isBigWig;          // flag indicates file is BigWig format;

    // BBFile Header items - Table C:
    // mMagic number (4 bytes) indicates file type and byte order :
    // 0x888FFC26 for bigWig, little endian if swapped
    // 0x8789F2EB for bigBed, little endian if swapped
    private int magic;                // 4 byte mMagic Number
    private short version;            // 2 byte version ID; currently 3
    private short nZoomLevels;         // 2 byte count of zoom sumary levels
    private long chromTreeOffset;     // 8 byte offset to mChromosome B+ Tree index
    private long fullDataOffset;      // 8 byte offset to unzoomed data dataCount
    private long fullIndexOffset;     // 8 byte offset to R+ Tree index of items
    private short fieldCount;         // 2 byte number of fields in bed. (0 for bigWig)
    private short definedFieldCount;  // 2 byte number of fields that are bed fields
    private long autoSqlOffset;       // 8 byte offset to 0 terminated string with .as spec
    private long totalSummaryOffset;  // 8 byte offset to file summary data block
    private int uncompressBuffSize;  // 4 byte maximum size for decompressed buffer
    private long reserved;            // 8 bytes reserved for future expansion. Currently 0

    // constructor reads BBFile header from an input stream
    public BBFileHeader(String path, SeekableStream fis, long fileOffset) {


        // save the path and seekable file handle
        this.path = path;
        this.fis = fis;
        fileHeaderOffset = fileOffset;

        // read in BBFile header
        isHeaderOK = readBBFileHeader(fileHeaderOffset);

    }

    /*
    *   Constructor loads BBFile header class from parameter specifications.
    *
    *   Parameters: (as defined above)
    * */
    public BBFileHeader(
            int magic,
            short version,
            short zoomLevels,
            long chromTreeOffset,
            long fullDataOffset,
            long fullIndexOffset,
            short fieldCount,
            short definedFieldCount,
            long autoSqlOffset,
            long totalSummaryOffset,
            int uncompressBuffSize,
            long reserved) {

        this.magic = magic;

        // Note: may want to validate the rest of the fields as well
        if (isBigWig() || isBigBed())
            this.isHeaderOK = true;

        this.version = version;
        this.nZoomLevels = zoomLevels;
        this.chromTreeOffset = chromTreeOffset;
        this.fullDataOffset = fullDataOffset;
        this.fullIndexOffset = fullIndexOffset;
        this.fieldCount = fieldCount;
        this.definedFieldCount = definedFieldCount;
        this.autoSqlOffset = autoSqlOffset;
        this.totalSummaryOffset = totalSummaryOffset;
        this.uncompressBuffSize = uncompressBuffSize;
        this.uncompressBuffSize = uncompressBuffSize;
        this.reserved = reserved;
    }

    public String getPath() {
        return path;
    }


    public boolean isHeaderOK() {
        return isHeaderOK;
    }


    public boolean isLowToHigh() {
        return isLowToHigh;
    }

    public boolean isBigBed() {
        return isBigBed;
    }

    public boolean isBigWig() {
        return isBigWig;
    }

    public int getFileHeaderSize() {
        return BBFILE_HEADER_SIZE;
    }

    // ************* return header items ****************

    public int getMagic() {
        return magic;
    }

    public short getVersion() {
        return version;
    }

    public short getZoomLevels() {
        return nZoomLevels;
    }

    public long getChromosomeTreeOffset() {
        return chromTreeOffset;
    }

    public long getFullDataOffset() {
        return fullDataOffset;
    }

    public long getFullIndexOffset() {
        return fullIndexOffset;
    }

     public short getFieldCount() {
         return fieldCount;
     }

     public short getDefinedFieldCount() {
         return definedFieldCount;
     }

    public long getAutoSqlOffset() {
        return autoSqlOffset;
    }

    public long getTotalSummaryOffset() {
        return totalSummaryOffset;
    }

    public int getUncompressBuffSize() {
        return uncompressBuffSize;
    }


    public void print() {

        if (isHeaderOK) {
            if (isBigWig())
                System.out.println("BigWig file " + path + ", file header at location " + fileHeaderOffset);
            else if (isBigBed())
                System.out.println("BigBed file " + path + ", file header at location " + fileHeaderOffset);
        } else {
            System.out.println("BBFile " + path + "  with bad magic = " + magic +
                    " from file header location " + fileHeaderOffset);
            return; // bad read - remaining header items not interpreted
        }


        // header fields
        System.out.println("BBFile header magic = " + magic);
        System.out.println("Version = " + version);
        System.out.println("Zoom Levels = " + nZoomLevels);
        System.out.println("Chromosome Info B+ tree offset = " + chromTreeOffset);
        System.out.println("Data Block offset = " + fullDataOffset);
        System.out.println("Chromosome Data R+ tree offset = " + fullIndexOffset);
        System.out.println("Bed fields count = " + fieldCount);
        System.out.println("Bed defined fields count = " + definedFieldCount);
        System.out.println("AutoSql Offset = " + autoSqlOffset);
        System.out.println("Total Summary offset = " + totalSummaryOffset);
        System.out.println("Maximum uncompressed buffer size = " + uncompressBuffSize);
        System.out.println("m_reserved = " + reserved);
    }

    /*
     *  Reads in BBFile header information.
     *
     *  Returns:
     *      Success status flag is true for successfully read header,
     *      or is false for a read error.
    **/
    private boolean readBBFileHeader(long fileOffset) {

        BBFileHeader bbHeader = null;
        LittleEndianInputStream lbdis = null;
        DataInputStream bdis = null;

        byte[] buffer = new byte[BBFILE_HEADER_SIZE];


        try {
            // Read bigBed header into a buffer
            fis.seek(fileOffset);
            fis.readFully(buffer);

            // decode header - determine byte order from first 4 bytes
            // first assume byte order is low to high
            isLowToHigh = true;
            lbdis = new LittleEndianInputStream(new ByteArrayInputStream(buffer));
            magic = lbdis.readInt();

            // check for a valid bigBed or bigWig file
            if (magic == BIGWIG_MAGIC_LTH)
                isBigWig = true;
            else if (magic == BIGBED_MAGIC_LTH)
                isBigBed = true;

                // try high to low byte order
            else {
                bdis = new DataInputStream(new ByteArrayInputStream(buffer));
                magic = bdis.readInt();

                // check for a valid bigBed or bigWig file
                if (magic == BIGWIG_MAGIC_HTL)
                    isBigWig = true;
                else if (magic == BIGBED_MAGIC_HTL)
                    isBigBed = true;

                else
                    return false;   // can't identify BBFile type

                // success - set order high to low
                isLowToHigh = false;
            }

            // Get header information
            if (isLowToHigh) {
                version = lbdis.readShort();
                nZoomLevels = lbdis.readShort();
                chromTreeOffset = lbdis.readLong();
                fullDataOffset = lbdis.readLong();
                fullIndexOffset = lbdis.readLong();
                fieldCount = lbdis.readShort();
                definedFieldCount = lbdis.readShort();
                autoSqlOffset = lbdis.readLong();
                totalSummaryOffset = lbdis.readLong();
                uncompressBuffSize = lbdis.readInt();
                reserved = lbdis.readLong();
            } else {
                version = bdis.readShort();
                nZoomLevels = bdis.readShort();
                chromTreeOffset = bdis.readLong();
                fullDataOffset = bdis.readLong();
                fullIndexOffset = bdis.readLong();
                fieldCount = bdis.readShort();
                definedFieldCount = bdis.readShort();
                autoSqlOffset = bdis.readLong();
                totalSummaryOffset = bdis.readLong();
                uncompressBuffSize = bdis.readInt();
                reserved = bdis.readLong();
            }

        } catch (IOException ex) {
            throw new RuntimeException("Error reading file header for " + path, ex);
        }

        // file header was read properly
        return true;
    }

}  // mEndBase of class BBFileHeader
