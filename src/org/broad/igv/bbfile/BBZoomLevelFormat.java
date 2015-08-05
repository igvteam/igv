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
 * Date: Jan 18, 2010
 * Time: 8:44:33 PM
 * To change this template use File | Settings | File Templates.
 */
/*
*   Container class for holding zoom level format information, BBFile Table O.
*
*   Note: The header section corresponds to Table O and the file location for
*   zoomCount is identified by the dataOffset entry in BBFile zoom level headers
*   Table D.
*
*   The zoomCount field in Table O is uncompressed, but the zoomData section is compressed.
*
*   Note:
*   Table P describes the zoom data record format. It contain statistics for
*   mChromosome regions in the zoom data. There are  Table O zoomCount data records
*   read into the ZoomDataRecords array.
*
* */
public class BBZoomLevelFormat {

    private static Logger log = Logger.getLogger(BBZoomLevelFormat.class);

    // size of zoom level format header is 4 byte; only accounts for zoomCount of Table O
    public static final int ZOOM_FORMAT_HEADER_SIZE = 4;

    // Perhaps there is a better way to know the top permitted record count
    public static final int MAX_ZOOM_DATA_RECORDS = 100000000;

    // defines the zoom data file access
    private int zoomLevel;             // zoom level for data
    private SeekableStream fis;      // BBFile handle
    private long zoomFormatOffset;     // BBFile zoom level data format offset
    private long zoomDataOffset;       // BBFile zoom level data offset
    private long zoomIndexOffset;      // BBFile zoom level R+ tree offset

    // data reader members
    private boolean isLowToHigh;   // zoom data low to high byte order if true; else high to low

    // zoom level data - BBFile Table O
    private int zoomRecordCount;   // number of data records; zoomCount in Table O
    private long zoomDataSize;     // number of (compressed) zoom level data bytes

    /*
    *   constructor   - reads zoom level format (but not the data) for a zoom level.
    *
    *   Parameters:
    *       zoomLevel - zoom level; from 1 to zoomLevels in Table C File Header
    *       fis - input stream reader handle
    *       fileOffset - file location for the zoom format table
    * `     datasize - byte size of (compressed/uncompressed) zoom data block
    *       isLowToHigh - boolean flag indicates if buffer data is arranged low to high byte
    *       uncompressBufSize - buffer size for decompressed data; or 0 for uncompressed
    *
    *   Note: Zoom level data Table O is arranged as:
    *       zoomCount - 4 bytes
    *       zoomData - dataSize bytes
    *       R+ zoom index starts at fileOffset + dataSize
    * */
    public BBZoomLevelFormat(int zoomLevel, SeekableStream fis, long fileOffset,
                           long dataSize, boolean isLowToHigh, int uncompressBufSize) {

        // store file access info
        this.zoomLevel = zoomLevel;
        this.fis = fis;
        zoomFormatOffset = fileOffset;
        zoomDataSize = dataSize;
        this.isLowToHigh = isLowToHigh;

        // Note: a bad zoom data header will result in a 0 count returned
        // or an IOException

        // size of buffer is ZOOM_FORMAT_HEADER_SIZE to get the record count
        byte[] buffer = new byte[ZOOM_FORMAT_HEADER_SIZE];

        try {

            // Read zoom level data format into a buffer
            fis.seek(zoomFormatOffset);
            fis.readFully(buffer);

            // decode header - or fail
            if(this.isLowToHigh) {
                LittleEndianInputStream lbdis = new LittleEndianInputStream(new ByteArrayInputStream(buffer));
                zoomRecordCount = lbdis.readInt();
            }
            else {
                DataInputStream bdis = new DataInputStream(new ByteArrayInputStream(buffer));
                zoomRecordCount = bdis.readInt();
            }

        } catch (IOException ex) {
            log.error("Error reading zoom level data records (Table O) ", ex);
            throw new RuntimeException("Error reading zoom level data records (Table O)", ex);
            }

        // integrity check - should be > 0 or less than a max like 100M records?
            // Note: if trouble reading zoom data records, readAllZoomLevelRecords returns 0
            if(zoomRecordCount < 0 || zoomRecordCount > MAX_ZOOM_DATA_RECORDS)
                return;  // terminate if bad zoom level data encountered

            // Position file offset past the current zoom level header to pick up
            // the zoom data records which immediately follow.
            zoomDataOffset = zoomFormatOffset + ZOOM_FORMAT_HEADER_SIZE;

            // calculate the position of the R+ zoom index tree
            zoomIndexOffset = zoomDataOffset + zoomDataSize;
    }

    /*
    *   Method returns the zoom level for zoom data.
    *
    *   Returns:
    *       zoom level for zoom data
    * */
    public int getZoomLevel() {
           return zoomLevel;
    }

    /*
    *   Method returns the file input stream handle.
    *
    *   Returns:
    *       file input stream handle
    * */
     public SeekableStream getBBFis() {
        return fis;
    }

    /*
    *   Method returns the file location for zoom data format.
    *
    *   Returns:
    *       file location for zoom data format
    * */
    public long getZoomFormatLocation() {
       return zoomFormatOffset;
    }

    /*
    *   Method returns the zoom count found in the zoom level format header .
    *
    *   Note: Zoom count is zoomCount in Table O.
    *
    *   Returns:
    *       file location for zoom level R+ data records
    * */
    public int getZoomRecordCount() {
        return zoomRecordCount;
    }

    /*
    *   Method returns the file location for zoom level data.
    *
    *   Note: ZoomFormatOffset (Table D dataOffset) + zoom data header size
    *
    *   Returns:
    *       file location for zoom level R+ data records
    * */
    public long getZoomDataOffset() {
        return zoomDataOffset;
    }

    /*
    *   Method returns the calculated compressed zoom level data byte size.
    *
    *   Note: Compressed data size is calculated from Table D
    *   (indexOffset - dataOffset) and is provided by the constructor.
    *
    *   Returns:
    *       calculated compressed zoom level data byte size
    * */
    public long getZoomDataSize() {
        return zoomDataSize;
    }

    /*
    *   Method returns the calculated file location for zoom level R+ tree.
    *
    *   Note: should be equal to indexOffset in zoom level header - Table D
    *
    *   Returns:
    *       file location for zoom level R+ tree
    * */
     public long getZoomIndexOffset() {
        return zoomIndexOffset;
    }


     public void print(){

        // identify the zoom level and record count
       log.debug("Zoom level " + zoomLevel + " format Table O found at "
             + zoomFormatOffset);
       log.debug("Zoom record count is " + zoomRecordCount);
       log.debug("Zoom data location is " + zoomDataOffset);
       log.debug("Zoom data size is " + zoomDataOffset);
       log.debug("Zoom index tree location is " + zoomIndexOffset);
    }

}
