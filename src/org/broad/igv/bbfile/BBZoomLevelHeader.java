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

import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: martind
 * Date: Dec 17, 2009
 * Time: 4:36:21 PM
 *
 * To change this template use File | Settings | File Templates.
 */
/*
*   Container class for holding zoom level header information, BBFile Table D.
*
*   Constructed either from BBFile read or by load of header values.
*
* */
public class BBZoomLevelHeader {

    private static Logger log = Logger.getLogger(BBZoomLevelHeader.class);

    static public final int ZOOM_LEVEL_HEADER_SIZE = 24;

    // Defines the Big Binary File (BBFile) access
    private SeekableStream fis;          // BBFile input stream handle
    private long zoomLevelHeaderOffset;    // file location for zoom level header
    int zoomLevel;     // the zoom level for this information

    // zoom level header information - BBFile Table D
    private int reductionLevel;   // number of bases summerized
    private int reserved;         // reserved, currently 0
    private long dataOffset;      // file position of zoom data
    private long indexOffset;     // file position for index of zoomed data

    /*
    *   Constructor reads zoom level header
    *
    *   Parameters:
    *       fis - File input stream handle
    *       fileOffset - file byte position for zoom header
    *       zoomLevel - level of zoom
    *       isLowToHigh - indicates byte order is low to high, else is high to low
    * */
    public BBZoomLevelHeader(SeekableStream fis, long fileOffset, int zoomLevel,
                             boolean isLowToHigh){

        this.fis = fis;
        zoomLevelHeaderOffset = fileOffset;
        this.zoomLevel = zoomLevel;

        readZoomLevelHeader(zoomLevelHeaderOffset, this.zoomLevel, isLowToHigh);
    }
    /*
    *   Constructor loads zoom level header according to parameter specification.
    *
    *   Parameters: (as defined above)
    * */
    public BBZoomLevelHeader(int zoomLevel, int reductionLevel, int reserved,
                             long dataOffset, long indexOffset){
        this.zoomLevel = zoomLevel;
        this.reductionLevel = reductionLevel;
        this.reserved = reserved;
        this.dataOffset = dataOffset;
        this.indexOffset = indexOffset;
    }

    /*
    *   Method returns the zoom level.
    *
    *   Returns:
    *       zoom level
    * */
    public int getZoomLevel() {
        return zoomLevel;
    }

    /*
    *   Method returns the reduction level for the zoom level.
    *
    *   Returns:
    *       reduction level
    * * */
    public int getReductionLevel() {
        return reductionLevel;
    }

    /*
    *   Method returns the reserved value.
    *
    *   Returns:
    *       reserved value
    * * */
    public int getReserved() {
        return reserved;
    }

    /*
    *   Method returns the zoom level data file location.
    *
    *   Returns:
    *       zoom level data file location
    * */
    public long getDataOffset() {
        return dataOffset;
    }

    /*
    *   Method returns the zoom level R+ index tree file location.
    *
    *   Returns:
    *       R+ index tree file location
    * */
    public long getIndexOffset() {
        return indexOffset;
    }

    /*
    *   Method prints the zoom level header info.
    * */
    public void print(){

        // Table D - Zoom Level Header information
        System.out.println("Zoom level " + zoomLevel + " header Table D: ");
        System.out.println("Number of zoom level bases = " + reductionLevel);
        System.out.println("Reserved = " + reserved);
        System.out.println("Zoom data offset = " + dataOffset);
        System.out.println("Zoom index offset = " + indexOffset);
    }

    /*
    *   Reads zoom level header information into class data members.
    *
    *   Parameters:
    *       fileOffset - Byte position in fle for zoom header
    *       zoomLevel - level of zoom
    *       isLowToHigh - indicate byte order is low to high, else is high to low
    * */
    private void readZoomLevelHeader(long fileOffset, int zoomLevel, boolean isLowToHigh) {

       LittleEndianInputStream lbdis = null;
       DataInputStream bdis = null;

        byte[] buffer = new byte[ZOOM_LEVEL_HEADER_SIZE];

            try {

            // Read zoom header into a buffer
            fis.seek(fileOffset);
            fis.readFully(buffer);

            // decode header
            if(isLowToHigh)
                lbdis = new LittleEndianInputStream(new ByteArrayInputStream(buffer));
            else
                bdis = new DataInputStream(new ByteArrayInputStream(buffer));

            // Get zoom level information
            if(isLowToHigh){
                reductionLevel = lbdis.readInt();
                reserved = lbdis.readInt();
                dataOffset = lbdis.readLong();
                indexOffset = lbdis.readLong();
            }
            else {
                reductionLevel = bdis.readInt();
                reserved = bdis.readInt();
                dataOffset = bdis.readLong();
                indexOffset = bdis.readLong();
            }

        }catch(IOException ex) {
            log.error("Error reading zoom level header: " + zoomLevel, ex);
            throw new RuntimeException("Error reading zoom header " + zoomLevel, ex);
        }
    }


}
