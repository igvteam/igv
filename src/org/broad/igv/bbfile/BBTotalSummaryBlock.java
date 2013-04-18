/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.bbfile;

import net.sf.samtools.seekablestream.SeekableStream;
import org.apache.log4j.Logger;
import org.broad.tribble.util.LittleEndianInputStream;

import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: martind
 * Date: Jan 18, 2010
 * Time: 12:51:25 PM
 * To change this template use File | Settings | File Templates.
 */
/*
*   Containeer class for file statistics - BBFile Table DD
*
* */
public class BBTotalSummaryBlock {

    private static Logger log = Logger.getLogger(BBTotalSummaryBlock.class);

    public static final int TOTAL_SUMMARY_BLOCK_SIZE = 40;

    // defines the R+ Tree access
    private SeekableStream fis;      // BBFile handle
    private long summaryBlockOffset;   // file offset to TotalSummaryBlock

    // File data statistics for calculating mean and standard deviation
    private long basesCovered;     // number of bases with data
    private float minVal;          // minimum value for file data
    private float maxVal;          // maximum value for file data
    private float sumData;         // sum of all squares of file data values
    private float sumSquares;      // sum of all squares of file data values

    /*
   *   Constructor for reading in TotalSummaryBlock from BBFile
   *
   *    Parameters:
   *    fis - file input stream handle
   *    fileOffset - file offset to TotalSummaryBlock
   *    isLowToHigh - indicates byte order is low to high if true, else is high to low
   * */
    public BBTotalSummaryBlock(SeekableStream fis, long fileOffset, boolean isLowToHigh)
    {

        LittleEndianInputStream lbdis = null;
        DataInputStream bdis = null;
        
        byte[] buffer = new byte[TOTAL_SUMMARY_BLOCK_SIZE];

        // save the seekable file handle  and B+ Tree file offset
        this.fis = fis;
        summaryBlockOffset = fileOffset;

        try {
            // Read TotalSummaryBlock header into a buffer
            fis.seek(fileOffset);
            fis.readFully(buffer);

            // decode header
            if(isLowToHigh)
                lbdis = new LittleEndianInputStream(new ByteArrayInputStream(buffer));
            else
                bdis = new DataInputStream(new ByteArrayInputStream(buffer));

            // Get TotalSummaryBlcok information
            if(isLowToHigh){
                basesCovered = lbdis.readLong();
                minVal = lbdis.readFloat();
                maxVal = lbdis.readFloat();
                sumData = lbdis.readFloat();
                sumSquares = lbdis.readFloat();
            }
            else {
                basesCovered = bdis.readLong();
                minVal = bdis.readFloat();
                maxVal = bdis.readFloat();
                sumData = bdis.readFloat();
                sumSquares = bdis.readFloat();
            }

        }catch(IOException ex) {
            log.error("Error reading Total Summary Block ", ex);
            throw new RuntimeException("Error reading Total Summary Block", ex);
            }

        }

    /*
    *   Constructor for filling in TotalSummaryBlock
    * */
    public BBTotalSummaryBlock(long basesCovered, float minVal, float maxVal,
                               float sumData, float sumSquares){

        this.basesCovered = basesCovered;
        this.minVal = minVal;
        this.maxVal = maxVal;
        this.sumData = sumData;
        this.sumSquares = sumSquares;

    }

    public static int getSummaryBlockSize() {
        return TOTAL_SUMMARY_BLOCK_SIZE;
    }

    public SeekableStream getMBBFis() {
        return fis;
    }

    public long getSummaryBlockOffset() {
        return summaryBlockOffset;
    }

     public long getBasesCovered() {
        return basesCovered;
    }

    public float getMinVal() {
        return minVal;
    }

    public float getMaxVal() {
        return maxVal;
    }

    public float getSumData() {
        return sumData;
    }

    public float getSumSquares() {
        return sumSquares;
    }

    public void printTotalSummaryBlock(){

        // Table D - Zoom Level Header information
        log.debug("BBFile TotalSummaryBlock (Table DD):");
        log.debug("Number of bases covered= " + basesCovered);
        log.debug("MinVal = " + minVal);
        log.debug("MaxVal = " + maxVal);
        log.debug("Sum of data values = "+ sumData);
        log.debug("Sum of squares values = " + sumSquares);
    }

}
