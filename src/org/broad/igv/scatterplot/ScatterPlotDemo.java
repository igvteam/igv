/* ---------------------
 * IGVScatterPlotDemo.java
 * ---------------------
 * (C) Copyright 2002-2009, by Object Refinery Limited.
 *
 */

package org.broad.igv.scatterplot;

//import demo.SampleXYDataset2;


import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
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
        igvPlotFrame = new ScatterPlotFrame("IGV Data: chr7: 83,584,980 - 123,460,500", igvData);

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
            cnData[i] = ((double) i / sampleCount) * 4 * (2*Math.random());
            exprData[i] = -2 + ((double) i / sampleCount) * 4 * (2*Math.random());
            methData[i] = 1 - i * Math.random() / sampleCount;
            System.out.println(cnData[i]+ "\t" + exprData[i] + "\t" + methData[i]);
        }

        Map<String, double[]> dataMap = new LinkedHashMap<String, double[]>();
        dataMap.put("Copy Number", cnData);
        dataMap.put("Expression", exprData);
        dataMap.put("Methylation", methData);

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
            System.out.println(cnData[i]+ "\t" + exprData[i] + "\t" + methData[i]);
       }

        Map<String, String[]> symbolMap = new LinkedHashMap<String, String[]>();
        symbolMap.put("Treatment", treated);
        symbolMap.put("Hyper mutated", hyperMutated);
        symbolMap.put("Cluster Type", cluster);

        // load the base class ScatterPlotData
        return new ScatterPlotData( samples,  symbolMap, dataMap);

    }

    /**
     * Starting point for the demonstration application.
     *
     * @param args  ignored.
     */
    public static void main(String[] args) {
        //int demoPoints = 100;

        int demoPoints = 100; //Integer.parseInt(args[0].trim());

        ScatterPlotDemo demo = new ScatterPlotDemo(demoPoints);

    }

}
