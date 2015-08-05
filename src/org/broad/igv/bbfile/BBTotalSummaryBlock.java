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
