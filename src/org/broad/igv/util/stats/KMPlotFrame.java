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

/*
 * Created by JFormDesigner on Sat Dec 04 19:26:01 EST 2010
 */

package org.broad.igv.util.stats;

import java.awt.*;
import java.awt.event.*;
import java.util.*;
import javax.swing.*;
import javax.swing.border.*;

import org.broad.igv.track.AttributeManager;
import org.broad.igv.track.Track;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * @author jrobinso
 */
public class KMPlotFrame extends JFrame {

    Collection<Track> tracks;
    private XYPlot plot;

    public KMPlotFrame(Collection<Track> tracks) {

        //setLocationRelativeTo(owner);

        this.tracks = tracks;

        initComponents();

        censureColumn.addItem("");
        sampleColumn.addItem("");
        survivalColumn.addItem("");
        groupByColumn.addItem("");
        for (String key : AttributeManager.getInstance().getAttributeKeys()) {
            censureColumn.addItem(key);
            sampleColumn.addItem(key);
            survivalColumn.addItem(key);
            groupByColumn.addItem(key);
        }

        censureColumn.setSelectedItem("CENSURED");
        sampleColumn.setSelectedItem("LINKING_ID");
        survivalColumn.setSelectedItem("SURVIVAL (DAYS)");
        groupByColumn.setSelectedItem("SUBTYPE");


        XYDataset dataset = updateDataset();
        JFreeChart chart = ChartFactory.createXYStepChart("Kaplan-Meier", "Days", "Survival", dataset,
                PlotOrientation.VERTICAL, true, true, false);
        plot = chart.getXYPlot();
        ChartPanel plotPanel = new ChartPanel(chart);

        contentPanel.add(plotPanel, BorderLayout.CENTER);

    }


    private void closeButtonActionPerformed(ActionEvent e) {
        setVisible(false);
    }

    public XYDataset updateDataset() {

        ArrayList<DataPoint> dataPoints = new ArrayList(tracks.size());
        HashSet<String> participants = new HashSet();
        for (Track t : tracks) {
            try {
                String part = t.getAttributeValue(sampleColumn.getSelectedItem().toString());
                if (!participants.contains(part)) {
                    participants.add(part);
                    String survivalString = t.getAttributeValue(survivalColumn.getSelectedItem().toString());
                    int survival = Integer.parseInt(survivalString);
                    String censureString = t.getAttributeValue(censureColumn.getSelectedItem().toString());
                    boolean censured = censureString != null && censureString.equals("1");
                    String group = t.getAttributeValue(groupByColumn.getSelectedItem().toString());
                    dataPoints.add(new DataPoint(part, survival, censured, group));
                } else {
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


        //XYSeries series1;
        XYSeriesCollection dataset = new XYSeriesCollection();
        for (Map.Entry<String, ArrayList<DataPoint>> entry : map.entrySet()) {

            java.util.List<DataPoint> pts = entry.getValue();
            Collections.sort(pts);

            int[] time = new int[pts.size()];
            boolean[] censured = new boolean[pts.size()];
            for (int i = 0; i < pts.size(); i++) {
                time[i] = Math.max(1, pts.get(i).time / 30);
                censured[i] = pts.get(i).censured;
            }

            java.util.List<KaplanMeierEstimator.Interval> controlIntervals = KaplanMeierEstimator.compute(time, censured);


            XYSeries series1 = new XYSeries(entry.getKey());
            for (KaplanMeierEstimator.Interval interval : controlIntervals) {
                series1.add(interval.getEnd(), interval.getCumulativeSurvival());
            }
            dataset.addSeries(series1);
        }

        return dataset;

    }

    private void updateButtonActionPerformed(ActionEvent e) {
        XYDataset dataset = updateDataset();
        plot.setDataset(dataset);
        repaint();

    }

    public static class DataPoint implements Comparable<DataPoint> {
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


    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        panel1 = new JPanel();
        censureColumn = new JComboBox();
        sampleColumn = new JComboBox();
        survivalColumn = new JComboBox();
        groupByColumn = new JComboBox();
        label1 = new JLabel();
        label2 = new JLabel();
        label3 = new JLabel();
        label4 = new JLabel();
        updateButton = new JButton();
        buttonBar = new JPanel();
        closeButton = new JButton();

        //======== this ========
        setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
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
                    panel1.setLayout(null);
                    panel1.add(censureColumn);
                    censureColumn.setBounds(145, 62, 215, censureColumn.getPreferredSize().height);
                    panel1.add(sampleColumn);
                    sampleColumn.setBounds(145, 0, 215, sampleColumn.getPreferredSize().height);
                    panel1.add(survivalColumn);
                    survivalColumn.setBounds(145, 31, 215, survivalColumn.getPreferredSize().height);
                    panel1.add(groupByColumn);
                    groupByColumn.setBounds(145, 93, 215, groupByColumn.getPreferredSize().height);

                    //---- label1 ----
                    label1.setText("Sample column");
                    panel1.add(label1);
                    label1.setBounds(5, 5, 115, label1.getPreferredSize().height);

                    //---- label2 ----
                    label2.setText("Survival column");
                    panel1.add(label2);
                    label2.setBounds(5, 36, 115, 16);

                    //---- label3 ----
                    label3.setText("Censur column");
                    panel1.add(label3);
                    label3.setBounds(5, 67, 115, 16);

                    //---- label4 ----
                    label4.setText("Group by");
                    panel1.add(label4);
                    label4.setBounds(5, 98, 115, 16);

                    //---- updateButton ----
                    updateButton.setText("Update Plot");
                    updateButton.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent e) {
                            updateButtonActionPerformed(e);
                        }
                    });
                    panel1.add(updateButton);
                    updateButton.setBounds(410, 92, 145, updateButton.getPreferredSize().height);

                    { // compute preferred size
                        Dimension preferredSize = new Dimension();
                        for (int i = 0; i < panel1.getComponentCount(); i++) {
                            Rectangle bounds = panel1.getComponent(i).getBounds();
                            preferredSize.width = Math.max(bounds.x + bounds.width, preferredSize.width);
                            preferredSize.height = Math.max(bounds.y + bounds.height, preferredSize.height);
                        }
                        Insets insets = panel1.getInsets();
                        preferredSize.width += insets.right;
                        preferredSize.height += insets.bottom;
                        panel1.setMinimumSize(preferredSize);
                        panel1.setPreferredSize(preferredSize);
                    }
                }
                contentPanel.add(panel1, BorderLayout.NORTH);
            }
            dialogPane.add(contentPanel, BorderLayout.CENTER);

            //======== buttonBar ========
            {
                buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                buttonBar.setLayout(new GridBagLayout());
                ((GridBagLayout) buttonBar.getLayout()).columnWidths = new int[]{0, 80};
                ((GridBagLayout) buttonBar.getLayout()).columnWeights = new double[]{1.0, 0.0};

                //---- closeButton ----
                closeButton.setText("Close");
                closeButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        closeButtonActionPerformed(e);
                    }
                });
                buttonBar.add(closeButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 0, 0), 0, 0));
            }
            dialogPane.add(buttonBar, BorderLayout.SOUTH);
        }
        contentPane.add(dialogPane, BorderLayout.CENTER);
        setSize(595, 510);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JPanel panel1;
    private JComboBox censureColumn;
    private JComboBox sampleColumn;
    private JComboBox survivalColumn;
    private JComboBox groupByColumn;
    private JLabel label1;
    private JLabel label2;
    private JLabel label3;
    private JLabel label4;
    private JButton updateButton;
    private JPanel buttonBar;
    private JButton closeButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables


}
