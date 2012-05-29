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

