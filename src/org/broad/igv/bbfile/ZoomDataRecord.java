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

/**
 * Created by IntelliJ IDEA.
 * User: martind
 * Date: Jan 18, 2010
 * Time: 2:24:44 PM
 * To change this template use File | Settings | File Templates.
 */
/*
*   Container class for holding zoom level statistics, BBFile Table P.
*
* */
public class ZoomDataRecord {

    private static Logger log = Logger.getLogger(ZoomDataRecord.class);

    public static final int RECORD_SIZE = 32;

    private int zoomLevel;         // zoom level associated with data
    private int recordNumber;      // record number

    // chromosome region statistics (useful for calculating mean and standard deviation)
    private String chromName;      // chromosome/contig name
    private int chromId;           // Numerical ID for mChromosome/contig
    private int chromStart;        // starting base position  (from 0)
    private int chromEnd;          // ending base position
    private int basesCovered;        // number of bases with data
    private float minVal;          // minimum value for file data
    private float maxVal;          // maximum value for file data
    private float sumData;         // sum of all squares of file data values
    private float sumSquares;      // sum of squares of file data values

    /*
    *   Constructor for filling in zoom data record class.
    *
    *   Parameters:
    *       zoomLevel - level of zoom
    *       recordNumber - record sequence number of multiple zoom level records
    *       chromName - chromosome/contig name
    *       chromId - mChromosome ID
    *       chromstart - starting base for zoom data region
    *       chromEnd - ending base for zoom data region
     *      validCount - number of bases in the region for which there is data
     *      minVal - minimum value in region
     *      maxVal - maximum value in region
     *      sumData - sum of all region data
     *      sumSquares - sum of the squares of all region data
     *
    * */
    public ZoomDataRecord(int zoomLevel, int recordNumber, String chromName, int chromId, int chromStart, int chromEnd,
            int validCount, float minVal, float maxVal, float sumData, float sumSquares ){

        this.zoomLevel = zoomLevel;
        this.recordNumber = recordNumber;
        this.chromName = chromName;
        this.chromId = chromId;
        this.chromStart = chromStart;
        this.chromEnd = chromEnd;
        this.basesCovered = validCount;
        this.minVal = minVal;
        this.maxVal = maxVal;
        this.sumData = sumData;
        this.sumSquares = sumSquares;
    }

    public int getZoomLevel() {
        return zoomLevel;
    }

    public int getRecordNumber() {
        return recordNumber;
    }

     public String getChromName() {
        return chromName;
    }

     public int getChromId() {
        return chromId;
    }

    public int getChromStart() {
        return chromStart;
    }

    public int getChromEnd() {
        return chromEnd;
    }

    public int getBasesCovered() {
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

    public float getMeanVal() {
        return basesCovered == 0 ? 0 : sumData / basesCovered;
    }

    public float getSumSquares() {
        return sumSquares;
    }

    public void print(){

        // Table P - zoom data record
       log.debug("Zoom data record (Table DD) number " + recordNumber +
               " for zoom level " + zoomLevel);
        log.debug("ChromName = " + chromName);
        log.debug("ChromId = " + chromId);
        log.debug("ChromStart = " + chromStart);
        log.debug("ChromEnd = " + chromEnd);
        log.debug("ValidCount = " + basesCovered);
        log.debug("MinVal = " + minVal);
        log.debug("MaxVal = " + maxVal);
        log.debug("Sum of data values = " + sumData);
        log.debug("Sum of squares values = " + sumSquares);
    }
}

