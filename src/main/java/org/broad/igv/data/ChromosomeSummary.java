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
