/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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

package org.broad.igv.sam;

import org.apache.commons.math.stat.StatUtils;
import org.broad.igv.sam.reader.AlignmentQueryReader;
import org.broad.igv.sam.reader.SamQueryReaderFactory;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;


/**
 * @author jrobinso
 * @date Jan 14, 2011
 */
public class PairedEndStats {

    private double averageInsertSize;
    private double medianInsertSize;
    private double insertSizeStdev;
    private static final int MAX_PAIRS = 10000;

    public static void main(String[] args) throws IOException {

        AlignmentQueryReader reader = SamQueryReaderFactory.getReader(args[0], false);
        PairedEndStats stats = compute(reader.iterator());
        reader.close();

        System.out.println(args[0] + "\t" + stats.averageInsertSize + "\t" + stats.medianInsertSize + "\t" + stats.insertSizeStdev);

    }


    public PairedEndStats(double averageInsertSize, double medianInsertSize, double insertSizeStdev) {
        this.averageInsertSize = averageInsertSize;
        this.medianInsertSize = medianInsertSize;
        this.insertSizeStdev = insertSizeStdev;
    }

    public static PairedEndStats compute(String bamFile) {
        AlignmentQueryReader reader = null;
        try {
            reader = SamQueryReaderFactory.getReader(bamFile, false);
            PairedEndStats stats = compute(reader.iterator());
            return stats;
        } catch (IOException e) {
            System.out.println("Error reading sam file");
            e.printStackTrace();
            return null;
        }
        finally {
            try {
                if (reader != null)
                    reader.close();
            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }


    }

    public static PairedEndStats compute(Iterator<Alignment> alignments) {


        double[] insertSizes = new double[MAX_PAIRS];
        int nPairs = 0;

        while (alignments.hasNext()) {


            Alignment al = alignments.next();

            if (isProperPair(al)) {
                insertSizes[nPairs] = Math.abs(al.getInferredInsertSize());
                nPairs++;
            }


            if (nPairs >= MAX_PAIRS) {
                break;
            }

        }


        double mean = StatUtils.mean(insertSizes, 0, nPairs);
        double median = StatUtils.percentile(insertSizes, 0, nPairs, 50);
        double stdDev = Math.sqrt(StatUtils.variance(insertSizes, 0, nPairs));

        PairedEndStats stats = new PairedEndStats(mean, median, stdDev);

        return stats;


    }


    static boolean isProperPair(Alignment alignment) {
        if (alignment.isMapped() && alignment.isPaired() && alignment.isProperPair() &&
                !alignment.isDuplicate() && alignment.getMappingQuality() > 0 && !alignment.isVendorFailedRead()) {
            ReadMate mate = alignment.getMate();
            boolean mateMapped = mate != null && mate.isMapped();
            boolean sameChromosome = mateMapped && mate.getChr().equals(alignment.getChr());
            return mateMapped && sameChromosome;
        }
        return false;
    }

    public double getAverageInsertSize() {
        return averageInsertSize;
    }

    public double getMedianInsertSize() {
        return medianInsertSize;
    }

    public double getInsertSizeStdev() {
        return insertSizeStdev;
    }
}

