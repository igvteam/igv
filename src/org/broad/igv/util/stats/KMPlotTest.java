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

package org.broad.igv.util.stats;

import org.broad.igv.track.AttributeManager;
import org.broad.igv.track.Track;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYStepRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import java.util.*;


public class KMPlotTest {


    public static void main(String[] args) {
        testPlot();

    }

    /**
     * Data from http://www.mas.ncl.ac.uk/~nmf16/teaching/mas3311/handout4.pdf
     */
    private static void testPlot() {
        int[] controlTime = {1, 1, 2, 2, 3, 4, 4, 5, 5, 8, 8, 8, 8, 11, 11, 12, 12};
        boolean[] controlCensured = new boolean[controlTime.length];
        Arrays.fill(controlCensured, false);
        int[] drugTime = {6, 6, 6, 6, 7, 9, 10, 11, 11, 16, 17, 19, 20, 25, 32, 32};
        boolean[] drugCensured = {false, false, false, true, false, true, true, false, true, false, true,
                true, true, false, true, true};

        List<KaplanMeierEstimator.Interval> controlIntervals = KaplanMeierEstimator.compute(controlTime, controlCensured);
        List<KaplanMeierEstimator.Interval> drugIntervals = KaplanMeierEstimator.compute(drugTime, drugCensured);

        XYSeries series1;
        XYSeriesCollection dataset = new XYSeriesCollection();

        series1 = new XYSeries("Control");
        for (KaplanMeierEstimator.Interval interval : controlIntervals) {
            series1.add(interval.getEnd(), interval.getCumulativeSurvival());
        }
        dataset.addSeries(series1);

        XYSeries series2 = new XYSeries("Drug");
        for (KaplanMeierEstimator.Interval interval : drugIntervals) {
            series2.add(interval.getEnd(), interval.getCumulativeSurvival());
        }
        dataset.addSeries(series2);

        JFreeChart plot = ChartFactory.createXYLineChart(
                "Kaplan-Meier Test",
                "Time",
                "Estiamte",
                dataset,
                PlotOrientation.VERTICAL,
                true,
                true,
                false);
        XYStepRenderer renderer = new XYStepRenderer();
        ((XYPlot) plot.getPlot()).setRenderer(renderer);

        ChartPanel plotPanel = new ChartPanel(plot);


        KMPlotFrame frame = new KMPlotFrame();
        frame.setPlotPanel(plotPanel);

        frame.setVisible(true);

    }

    public static void openPlot(Collection<Track> tracks) {

        String participantColumn = "Linking_id";
        String survivalColumn = "Survival (days)";
        String censureColumn = "Censured";
        String groupByColumn = "Subtype";

        ArrayList<DataPoint> dataPoints = new ArrayList(tracks.size());
        HashSet<String> participants = new HashSet();
        for (Track t : tracks) {
            try {
                String part = t.getAttributeValue(participantColumn.toUpperCase());
                if (!participants.contains(part)) {
                    participants.add(part);
                    String survivalString = t.getAttributeValue(survivalColumn.toUpperCase());
                    int survival = Integer.parseInt(survivalString);
                    String censureString = t.getAttributeValue(censureColumn.toUpperCase());
                    boolean censured = censureString != null && censureString.equals("1");
                    String group = t.getAttributeValue(groupByColumn.toUpperCase());
                    dataPoints.add(new DataPoint(part, survival, censured, group));
                }
                else {
                    // TODO -- check consistency of participant data
                }
            } catch (NumberFormatException e) {
                // Just skip
            }
        }


        // Segregate by group
        Map<String, ArrayList<DataPoint>> map = new HashMap();
        for (DataPoint dp : dataPoints) {
            String g = dp.group;
            ArrayList<DataPoint> pts = map.get(g);
            if (pts == null) {
                pts = new ArrayList();
                map.put(g, pts);
            }
            pts.add(dp);
        }


        int maxEnd = 50;
        XYSeries series1;
        XYSeriesCollection dataset = new XYSeriesCollection();
        for (Map.Entry<String, ArrayList<DataPoint>> entry : map.entrySet()) {

            List<DataPoint> pts = entry.getValue();
            Collections.sort(pts);

            int[] time = new int[pts.size()];
            boolean[] censured = new boolean[pts.size()];
            for (int i = 0; i < pts.size(); i++) {
                time[i] = Math.max(1, pts.get(i).time / 30);
                censured[i] = pts.get(i).censured;
            }

            List<KaplanMeierEstimator.Interval> controlIntervals = KaplanMeierEstimator.compute(time, censured);


            series1 = new XYSeries(entry.getKey());
            for (KaplanMeierEstimator.Interval interval : controlIntervals) {
                series1.add(interval.getEnd(), interval.getCumulativeSurvival());
            }
            dataset.addSeries(series1);
        }


        JFreeChart plot = ChartFactory.createXYLineChart("Kaplan-Meier", "Months", "Survival", dataset,
                PlotOrientation.VERTICAL, true, true, false);
        XYStepRenderer renderer = new XYStepRenderer();
        ((XYPlot) plot.getPlot()).setRenderer(renderer);

        ChartPanel plotPanel = new ChartPanel(plot);


        KMPlotFrame frame = new KMPlotFrame();
        frame.setPlotPanel(plotPanel);

        frame.setVisible(true);


    }

    static class DataPoint implements Comparable<DataPoint> {
        String participant;
        int time;
        boolean censured;
        String group;

        DataPoint(String participant, int time, boolean censured, String group) {
            this.censured = censured;
            this.participant = participant;
            this.group = group;
            this.time = time;
        }

        public int compareTo(DataPoint dataPoint) {
            return time - dataPoint.time;
        }
    }

}
