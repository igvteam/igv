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

package org.broad.igv.charts;


import org.broad.igv.track.TrackType;
import org.broad.igv.ui.util.UIUtilities;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * A demo scatter plot.
 */
public class ScatterPlotDemo /*extends ApplicationFrame*/ {

    /**
     * A demonstration application showing a scatter plot.
     *
     * @param title  the frame title.
     */

    private static int igvSampleCount;              // count of IGV data points
    private static ScatterPlotData igvData;      // container class for IGV Data
    private static ScatterPlotFrame igvPlotFrame;     // scatterplot frame

    public ScatterPlotDemo(int sampleCount) {

        // create test data
        igvData = getTestData(sampleCount);

         // create scatteplot and plot test data
        igvPlotFrame = new ScatterPlotFrame(igvData);

        UIUtilities.invokeOnEventThread(new Runnable() {
            public void run() {
                igvPlotFrame.setVisible(true);
            }
        });
    }


    public ScatterPlotData getTestData(int sampleCount) {


        // sample data types (numeric measurements)
        String [] samples = new String[sampleCount];    // test sample names
        double [] cnData = new double[sampleCount];
        double [] exprData = new double[sampleCount];
        double [] methData = new double[sampleCount];

        for(int i=0; i<sampleCount; i++) {
            samples[i] = "Sample " + i;
//            cnData[i] = ((double) i / sampleCount) * 4 * (2*Math.random());
//            exprData[i] = -2 + ((double) i / sampleCount) * 4 * (2*Math.random());
//            methData[i] = 1 - i * Math.random() / sampleCount;
            double x = ((i - sampleCount/2) * 3.1415) / sampleCount;
            double y = 5* Math.sin(i);
            cnData[i] = x;
            exprData[i] = y;
            methData[i] = y;
        }

        Map<String, double[]> dataMap = new LinkedHashMap<String, double[]>();
        dataMap.put(TrackType.COPY_NUMBER.toString(), cnData);
        dataMap.put(TrackType.GENE_EXPRESSION.toString(), exprData);
        dataMap.put(TrackType.DNA_METHYLATION.toString(), methData);

        // sample attributes (categorical)
        String [] treated = new String[sampleCount];
        String [] hyperMutated = new String[sampleCount];
        String [] cluster = new String[sampleCount];

        for(int i=0; i<sampleCount; i++) {
            treated[i] = Math.random() > 0.5 ? "AA" : "B";
            if(i%6 == 0)
             treated[i] = "null";
            hyperMutated[i] = Math.random() > 0.9 ? "Y" : "N";
            if(i%7 == 0)
             hyperMutated[i] = "null";

            double classValue =  Math.random();
            if(classValue > .9)
                cluster[i] = "Proneural";
            else if(classValue > .8)
                cluster[i] = "neural";
            else if(classValue > .6)
                cluster[i] = "classical";
            else if(classValue > .4)
                cluster[i] = "mesenchymal";
            else
                cluster[i] = "unknown";
       }

        Map<String, String[]> symbolMap = new LinkedHashMap<String, String[]>();
        symbolMap.put("Treatment", treated);
        symbolMap.put("Hyper mutated", hyperMutated);
        symbolMap.put("Cluster Type", cluster);

        int [] mutationCount = new int[sampleCount];
         for(int i=0; i<sampleCount; i++) {
             mutationCount[i] = (int) (2 * Math.random());
         }

        // load the base class ScatterPlotData
        return new ScatterPlotData("Test", samples,  symbolMap, dataMap, mutationCount);

    }

    /**
     * Starting point for the demonstration application.
     *
     * @param args  ignored.
     */
    public static void main(String[] args) {
        //int demoPoints = 100;

        int demoPoints = 500; //Integer.parseInt(args[0].trim());

        ScatterPlotDemo demo = new ScatterPlotDemo(demoPoints);

    }

}
