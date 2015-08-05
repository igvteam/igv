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
 * Created by JFormDesigner on Sat Dec 04 19:26:01 EST 2010
 */

package org.broad.igv.util.stats;

import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.util.List;
import javax.swing.*;
import javax.swing.border.*;

import org.broad.igv.track.AttributeManager;
import org.broad.igv.track.Track;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYStepRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * @author jrobinso
 */
public class KMPlotFrame extends JFrame {

    Collection<Track> tracks;
    private XYPlot plot;

    // TODO -- control to set this
    int maxTime = 60;

    public KMPlotFrame(Collection<Track> tracks) {

        //setLocationRelativeTo(owner);

        this.tracks = tracks;

        initComponents();


        XYDataset dataset = updateDataset();
        JFreeChart chart = ChartFactory.createXYLineChart(
                "",
                "Time",
                "Survival",
                dataset,
                PlotOrientation.VERTICAL,
                true,
                true,
                false);

        XYStepRenderer renderer = new XYStepRenderer();
        plot = chart.getXYPlot();
        plot.setRenderer(renderer);
        ChartPanel plotPanel = new ChartPanel(chart);
        contentPanel.add(plotPanel, BorderLayout.CENTER);


        censurColumnControl.addItem("");
        //sampleColumnControl.addItem("");
        survivalColumnControl.addItem("");
        groupByControl.addItem("");

        final List<String> allAttributes = AttributeManager.getInstance().getAttributeNames();
        final List<String> groupableAttributes = AttributeManager.getInstance().getGroupableAttributes();

        // Populate pulldowns, and make guesses for column names.
        String survivalColumn = null;
        String censureColumn = null;
        String sampleColumn = null;
        for (String key : allAttributes) {
            censurColumnControl.addItem(key);
            //sampleColumnControl.addItem(key);
            survivalColumnControl.addItem(key);

            String tmp = key.toLowerCase();
            if (tmp.contains("survival") || tmp.contains("daystodeath")) survivalColumn = key;
            if (tmp.contains("censure")) censureColumn = key;
            if (tmp.contains("sample")) sampleColumn = key;
        }
        for (String key : groupableAttributes) {
            groupByControl.addItem(key);
        }


        if (survivalColumn != null) {
            survivalColumnControl.setSelectedItem(survivalColumn);
        }
        if (censureColumn != null) {
            censurColumnControl.setSelectedItem(censureColumn);
        }
        if (sampleColumn != null) {
            //sampleColumnControl.setSelectedItem(sampleColumn);
        }


    }


    private void closeButtonActionPerformed(ActionEvent e) {
        setVisible(false);
    }

    public XYDataset updateDataset() {

        XYSeriesCollection dataset = new XYSeriesCollection();

        final String survivalColumn = (String) survivalColumnControl.getSelectedItem();
        final String censureColumn = (String) censurColumnControl.getSelectedItem();
        final String groupByColumn = (String) groupByControl.getSelectedItem();
        if (survivalColumn != null) {
            ArrayList<DataPoint> dataPoints = new ArrayList(tracks.size());
            HashSet<String> participants = new HashSet();
            for (Track t : tracks) {
                try {
                    // Get the participant (sample) attribute value for this track
                    //final Object selectedItem = sampleColumnControl.getSelectedItem();
                    //if (selectedItem != null) {
                    String participant = t.getSample(); // t.getAttributeValue(selectedItem.toString());


                    if (!participants.contains(participant)) {   // Don't add same participant twice.
                        participants.add(participant);

                        // Get the survival time.
                        String survivalString = t.getAttributeValue(survivalColumn);
                        int survivalDays = Integer.parseInt(survivalString);
                        int survival = survivalDays;

                        // Is the patient censured at the end of the survival period?
                        String censureString = censureColumn == null ? null : t.getAttributeValue(censureColumn);
                        boolean censured = censureString != null && censureString.equals("1");

                        String group = groupByColumn == null ? null : t.getAttributeValue(groupByColumn);
                        if (group == null) group = "<No value>";
                        dataPoints.add(new DataPoint(participant, survival, censured, group));
                    } else {
                        // TODO -- check consistency of participant data
                    }
                    // }
                } catch (NumberFormatException e) {
                    // Just skip
                }
            }

            // Segregate by group
            Map<String, ArrayList<DataPoint>> map = new HashMap();
            for (DataPoint dp : dataPoints) {
                String g = dp.getGroup();
                ArrayList<DataPoint> pts = map.get(g);
                if (pts == null) {
                    pts = new ArrayList();
                    map.put(g, pts);
                }
                pts.add(dp);
            }


            //XYSeries series1;
            for (Map.Entry<String, ArrayList<DataPoint>> entry : map.entrySet()) {

                java.util.List<DataPoint> pts = entry.getValue();
                Collections.sort(pts);

                int[] time = new int[pts.size()];
                boolean[] censured = new boolean[pts.size()];
                for (int i = 0; i < pts.size(); i++) {
                    //int months = Math.max(1, pts.get(i).time / 30);  // <=  TODO -- HARDCODED MONTH DATE
                    time[i] = pts.get(i).time;
                    censured[i] = pts.get(i).censured;
                }

                java.util.List<KaplanMeierEstimator.Interval> controlIntervals = KaplanMeierEstimator.compute(time, censured);

                XYSeries series1 = new XYSeries(entry.getKey());
                for (KaplanMeierEstimator.Interval interval : controlIntervals) {
                    series1.add(interval.getEnd(), interval.getCumulativeSurvival());
                }
                dataset.addSeries(series1);
            }
        }
        return dataset;

    }

    private void survivalColumnControlActionPerformed(ActionEvent e) {
        XYDataset dataset = updateDataset();
        plot.setDataset(dataset);
        repaint();

    }

    public static class DataPoint implements Comparable<DataPoint> {
        String participant;
        int time;
        boolean censured;
        private String group;

        DataPoint(String participant, int time, boolean censured, String group) {
            this.censured = censured;
            this.participant = participant;
            this.group = group;
            this.time = time;
        }

        public int compareTo(DataPoint dataPoint) {
            return time - dataPoint.time;
        }

        public String getGroup() {
            return group;
        }

    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        panel1 = new JPanel();
        panel2 = new JPanel();
        label2 = new JLabel();
        survivalColumnControl = new JComboBox();
        panel3 = new JPanel();
        label3 = new JLabel();
        censurColumnControl = new JComboBox();
        panel4 = new JPanel();
        label4 = new JLabel();
        groupByControl = new JComboBox();

        //======== this ========
        setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("Kaplan-Meier Plot");
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(new BorderLayout());

                //======== panel1 ========
                {
                    panel1.setAlignmentX(0.0F);
                    panel1.setLayout(new BoxLayout(panel1, BoxLayout.Y_AXIS));

                    //======== panel2 ========
                    {
                        panel2.setAlignmentX(1.0F);
                        panel2.setLayout(null);

                        //---- label2 ----
                        label2.setText("Survival column");
                        panel2.add(label2);
                        label2.setBounds(new Rectangle(new Point(5, 10), label2.getPreferredSize()));

                        //---- survivalColumnControl ----
                        survivalColumnControl.addActionListener(new ActionListener() {
                            @Override
                            public void actionPerformed(ActionEvent e) {
                                survivalColumnControlActionPerformed(e);
                            }
                        });
                        panel2.add(survivalColumnControl);
                        survivalColumnControl.setBounds(120, 5, 235, survivalColumnControl.getPreferredSize().height);

                        { // compute preferred size
                            Dimension preferredSize = new Dimension();
                            for(int i = 0; i < panel2.getComponentCount(); i++) {
                                Rectangle bounds = panel2.getComponent(i).getBounds();
                                preferredSize.width = Math.max(bounds.x + bounds.width, preferredSize.width);
                                preferredSize.height = Math.max(bounds.y + bounds.height, preferredSize.height);
                            }
                            Insets insets = panel2.getInsets();
                            preferredSize.width += insets.right;
                            preferredSize.height += insets.bottom;
                            panel2.setMinimumSize(preferredSize);
                            panel2.setPreferredSize(preferredSize);
                        }
                    }
                    panel1.add(panel2);

                    //======== panel3 ========
                    {
                        panel3.setAlignmentX(1.0F);
                        panel3.setLayout(null);

                        //---- label3 ----
                        label3.setText("Censored column");
                        panel3.add(label3);
                        label3.setBounds(new Rectangle(new Point(5, 10), label3.getPreferredSize()));

                        //---- censurColumnControl ----
                        censurColumnControl.addActionListener(new ActionListener() {
                            @Override
                            public void actionPerformed(ActionEvent e) {
                                survivalColumnControlActionPerformed(e);
                            }
                        });
                        panel3.add(censurColumnControl);
                        censurColumnControl.setBounds(120, 5, 235, censurColumnControl.getPreferredSize().height);

                        { // compute preferred size
                            Dimension preferredSize = new Dimension();
                            for(int i = 0; i < panel3.getComponentCount(); i++) {
                                Rectangle bounds = panel3.getComponent(i).getBounds();
                                preferredSize.width = Math.max(bounds.x + bounds.width, preferredSize.width);
                                preferredSize.height = Math.max(bounds.y + bounds.height, preferredSize.height);
                            }
                            Insets insets = panel3.getInsets();
                            preferredSize.width += insets.right;
                            preferredSize.height += insets.bottom;
                            panel3.setMinimumSize(preferredSize);
                            panel3.setPreferredSize(preferredSize);
                        }
                    }
                    panel1.add(panel3);

                    //======== panel4 ========
                    {
                        panel4.setAlignmentX(1.0F);
                        panel4.setLayout(null);

                        //---- label4 ----
                        label4.setText("Group by");
                        panel4.add(label4);
                        label4.setBounds(new Rectangle(new Point(5, 10), label4.getPreferredSize()));

                        //---- groupByControl ----
                        groupByControl.addActionListener(new ActionListener() {
                            @Override
                            public void actionPerformed(ActionEvent e) {
                                survivalColumnControlActionPerformed(e);
                            }
                        });
                        panel4.add(groupByControl);
                        groupByControl.setBounds(120, 5, 235, groupByControl.getPreferredSize().height);

                        { // compute preferred size
                            Dimension preferredSize = new Dimension();
                            for(int i = 0; i < panel4.getComponentCount(); i++) {
                                Rectangle bounds = panel4.getComponent(i).getBounds();
                                preferredSize.width = Math.max(bounds.x + bounds.width, preferredSize.width);
                                preferredSize.height = Math.max(bounds.y + bounds.height, preferredSize.height);
                            }
                            Insets insets = panel4.getInsets();
                            preferredSize.width += insets.right;
                            preferredSize.height += insets.bottom;
                            panel4.setMinimumSize(preferredSize);
                            panel4.setPreferredSize(preferredSize);
                        }
                    }
                    panel1.add(panel4);
                }
                contentPanel.add(panel1, BorderLayout.NORTH);
            }
            dialogPane.add(contentPanel, BorderLayout.CENTER);
        }
        contentPane.add(dialogPane, BorderLayout.CENTER);
        setSize(565, 510);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JPanel panel1;
    private JPanel panel2;
    private JLabel label2;
    private JComboBox survivalColumnControl;
    private JPanel panel3;
    private JLabel label3;
    private JComboBox censurColumnControl;
    private JPanel panel4;
    private JLabel label4;
    private JComboBox groupByControl;
    // JFormDesigner - End of variables declaration  //GEN-END:variables


}
