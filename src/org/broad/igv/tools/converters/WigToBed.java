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

package org.broad.igv.tools.converters;

import org.broad.igv.tools.parsers.DataConsumer;
import org.broad.igv.tools.parsers.WiggleParser;
import org.broad.igv.track.TrackType;

import java.io.*;

/**
 * Converts a "wig" file to a "bed" file by applying a threshold
 *
 * @author jrobinso
 * @date Jan 20, 2011
 */
public class WigToBed implements DataConsumer {

    PrintWriter bedWriter;
    double lowerThreshold = 0.1;
    double higherThreshold = 0.3;
    String chr = null;
    int featureStart = -1;
    int featureEnd = -1;
    String type;
    private float score;

    static String upperName = "upper";
    static String lowerName = "lower";


    public static void main(String[] args) throws IOException {
        String input =  args[0]; //"/Users/jrobinso/Sigma/566.wgs.bam.large_isize.wig";
        WigToBed wigToBed = new WigToBed(input, .17f, .55f);
        WiggleParser parser = new WiggleParser(input, wigToBed, null);
        parser.parse();
    }

    public static void run(String inputFile, float hetThreshold, float homThreshold) throws IOException {

        lowerName = String.valueOf((int) (hetThreshold * 10));
        upperName =  String.valueOf((int) (homThreshold * 10));

        int len = inputFile.length();
        String outputFile = inputFile.substring(0, len - 4) + ".bed";
        WigToBed wigToBed = new WigToBed(outputFile, hetThreshold, homThreshold);
        WiggleParser parser = new WiggleParser(inputFile, wigToBed, null);
        parser.parse();
    }

    public WigToBed(String bedFile, float lowerThreshold, float higherThreshold) {
        this.lowerThreshold = lowerThreshold;
        this.higherThreshold = higherThreshold;
        try {
            bedWriter = new PrintWriter(new BufferedWriter(new FileWriter(bedFile)));
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }


    public void addData(String chr, int start, int end, float[] data, String name) {

        if (featureStart >= 0) {
            // Feature in progress
            if (start > featureEnd || data[0] < lowerThreshold) {
                writeCurrentFeature();
                featureStart = -1;
                score = 0f;
            }
            featureEnd = end;
        }

        if (featureStart >= 0) {

            if (data[0] > higherThreshold) {
                type =  String.valueOf(higherThreshold);
                score = Math.max(score, data[0]);

            } else if (data[0] < lowerThreshold) {
                writeCurrentFeature();
                featureStart = -1;
                score = 0f;
            }

        } else {
            if (data[0] > lowerThreshold) {
                featureStart = start;
                featureEnd = end;
                this.chr = chr;
                score = Math.max(score, data[0]);
                type = data[0] > higherThreshold ? String.valueOf(higherThreshold) : String.valueOf(lowerThreshold);
            }

        }
    }

    private void writeCurrentFeature() {
        bedWriter.println(this.chr + "\t" + featureStart + "\t" + featureEnd + "\t" + type);
        featureStart = -1;
    }

    public void parsingComplete() {
        if (featureStart >= 0) writeCurrentFeature();
        bedWriter.close();
    }

    public void setTrackParameters(TrackType trackType, String trackLine, String[] trackNames) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void setTrackParameters(TrackType trackType, String trackLine, String[] trackNames, boolean b) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * Set a tolerance for "sortedness" of the data.  A start position can be less than
     * the immediately previous start position by this amount.  This is needed for
     * chip-seq processing where the start position of an alignment can be artificially
     * extended post sorting.
     */
    public void setSortTolerance(int tolerance) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void setAttribute(String key, String value) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void setType(String type) {
        //To change body of implemented methods use File | Settings | File Templates.
    }


}
