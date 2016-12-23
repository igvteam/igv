/*
 *  The MIT License (MIT)
 *
 * Copyright (c) 2016 University of California San Diego
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 *
 */

package org.broad.igv.sam;

import org.apache.commons.math.stat.StatUtils;

/**
 * Created by jrobinso on 12/20/16.
 * <p>
 * Some simple statistics to distinguish dna, rna, and 3rd gen alignment datasets.
 */


public class ReadStats {

    int binSize = 1000;
    int nReads;
    int nOperatorCount = 0;
    int indelCount = 0;
    double medianReadLength = 0;
    double medianRefToReadRatio = 0;

    int nWorkingReads = 0;
    int nWorkingOperatorCount = 0;
    double[] readLengths;
    double[] refToReadRatios;

    public ReadStats(int binSize) {
        this.binSize = binSize;
        this.nReads = 0;
        this.nWorkingReads = 0;
        this.nWorkingOperatorCount = 0;
        this.readLengths = new double[binSize];
        this.refToReadRatios = new double[binSize];
    }

    public void addAlignment(Alignment alignment) {

        if (alignment.isMapped()) {

            final int readLength = alignment.getReadSequence().length();
            final char [] cigarString = alignment.getCigarString().toCharArray();

            readLengths[nWorkingReads] = readLength;
            nWorkingReads++;

            if (containsChar(cigarString, 'N')) {
                refToReadRatios[nWorkingOperatorCount] = ((double) (alignment.getAlignmentEnd() - alignment.getAlignmentStart())) / readLength;
                nWorkingOperatorCount++;
            }

            if (containsChar(cigarString, 'D') || containsChar(cigarString, 'I')) {
                indelCount++;
            }

            if (nWorkingReads == binSize) {
                compute();
            }
        }

    }

    public void compute() {

        if (nWorkingReads > 0) {
            double rl = StatUtils.percentile(readLengths, 0, nWorkingReads, 50);
            double f = ((double) nReads) / (nReads + nWorkingReads);
            medianReadLength = (f * medianReadLength + (1 - f) * rl);
        }

        if (nWorkingOperatorCount > 0) {
            double rr = StatUtils.percentile(refToReadRatios, 0, nWorkingOperatorCount, 50);
            double f = ((double) nOperatorCount) / (nOperatorCount + nWorkingOperatorCount);
            medianRefToReadRatio = (f * medianRefToReadRatio + (1 - f) * rr);
        }

        nReads += nWorkingReads;
        nOperatorCount += nWorkingOperatorCount;

        nWorkingReads = 0;
        nWorkingOperatorCount = 0;
        readLengths = new double[binSize];
        refToReadRatios = new double[binSize];

//        System.out.println("Median read length= " + medianReadLength);
//        System.out.println("Median (ref length / read length) ratio = " + medianRefToReadRatio);
//        System.out.println("# of reads with \"N\" operators : " + nOperatorCount +  "  (" + ((100.0 * nOperatorCount) / nReads) + "%");
//        System.out.println("# of reads with indels : " + indelCount +  "  (" + ((100.0 * indelCount) / nReads) + "%");

    }

    // Hand coded for speed
    private boolean containsChar(char [] charArray, char c) {
        for(int i=0; i<charArray.length; i++) {
            if(charArray[i] == c) return true;
        }
        return false;
    }


}
