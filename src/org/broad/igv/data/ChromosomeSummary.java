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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.data;

/**
 * @author jrobinso
 */
public class ChromosomeSummary {

    /**
     * The chromosome name
     */
    private String name;
    /**
     * The approximate number of data points on this chromsome.  Used to estimate data
     * loading completion %.
     */
    private int nDataPoints;
    /**
     * The approximate number of characters to the starting position for this
     * chromosomes data in the copy number file.
     */
    private long startPosition;

    /**
     * Creates a new instance of ChromsomeSummary
     */
    public ChromosomeSummary(String name, long startPosition) {

        // Convert name to ucsc convention
        this.name = name; //.startsWith("chr") ? name : "chr" + name;
        this.startPosition = startPosition;
    }

    public String getName() {
        return name;
    }

    public long getStartPosition() {
        return startPosition;
    }

    public int getNDataPts() {
        return nDataPoints;
    }

    public void setNDataPoints(int nSnps) {
        this.nDataPoints = nSnps;
    }


}
