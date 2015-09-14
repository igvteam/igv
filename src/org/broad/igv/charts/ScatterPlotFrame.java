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
 * Created by JFormDesigner on Wed Oct 26 22:35:11 EDT 2011
 */

package org.broad.igv.charts;

import org.broad.igv.track.TrackType;

import java.awt.*;
import java.awt.event.*;
import java.text.DecimalFormat;
import java.util.List;
import javax.swing.*;
import javax.swing.border.*;

/**
 * @author Stan Diamond
 */
public class ScatterPlotFrame extends JFrame {

    static String lastClassifySelection = null;

    ScatterPlotData scatterPlotData;
    private boolean deferUpdate;

    public ScatterPlotFrame() {
        initComponents();
    }

    public ScatterPlotFrame(ScatterPlotData scatterPlotData) {

        deferUpdate = true;

        this.scatterPlotData = scatterPlotData;

        setTitle(scatterPlotData.getTitle());

        initComponents();

        List<String> categoryList = scatterPlotData.getCategories();

        // Add data names (e.g. copy #, expression)
        for (String dn : scatterPlotData.getDataNames()) {
            categoryList.add(dn);
        }

        if (categoryList.size() > 0) {
            String[] categories = categoryList.toArray(new String[categoryList.size()]);
            classifyComboBox.setModel(new DefaultComboBoxModel(categories));
        } else {
            classifyComboBox.setEnabled(false);
        }


        List<String> dataTypeList = scatterPlotData.getDataNames();
        String[] dataTypes = dataTypeList.toArray(new String[dataTypeList.size()]);
        xAxisComboBox.setModel(new DefaultComboBoxModel(dataTypes));
        yAxisComboBox.setModel(new DefaultComboBoxModel(dataTypes));

        if (classifyComboBox.getItemCount() > 0) {
            if (lastClassifySelection != null) {
                classifyComboBox.setSelectedItem(lastClassifySelection);
            } else {
                classifyComboBox.setSelectedIndex(0);
            }
        }

        xAxisComboBox.setSelectedIndex(0);
        yAxisComboBox.setSelectedIndex(dataTypes.length > 1 ? 1 : 0);

        deferUpdate = false;
        updateModel();
    }


    private void closeMenuItemActionPerformed(ActionEvent e) {
        setVisible(false);
    }

    private void axisChanged(ActionEvent e) {
        updateModel();
    }

    private void attributeChanged(ActionEvent e) {
        lastClassifySelection = (String) classifyComboBox.getSelectedItem();
        updateModel();
    }

    /*
    *   createChart - method creates a JFree ChartPanel (derived from JPanel)
    *   which plots the XYDataset assigned form ScatterPlotData selections.
    *
    *   Parameters:
    *       title - string contains title of plot
    *       igvData - ScatterPlotFrame container class for IGV data samples
    *       xAxisName - IGV data measurement name for x Axis coordinate;
    *           must be a key defined in igvData dataMap.
    *       yAxisName - IGV data measurement name for y Axis coordinate;
    *           must be a key defined in igvData dataMap.
    *       attribute - sample attribute selected for series symbol display
    *       symbolSettings - symbol settings for series symbol display.
    *
    *   Returns:
    *       Scatterplot chart
    *
    * */


    private void updateModel() {

        if (deferUpdate) return;

        String xAxisName = (String) xAxisComboBox.getSelectedItem();
        String yAxisName = (String) yAxisComboBox.getSelectedItem();

        String[] sampleNames = scatterPlotData.getSampleNames();
        double[] xValues = scatterPlotData.getDataValues(xAxisName);
        double[] yValues = scatterPlotData.getDataValues(yAxisName);
        int[] mutationCount = scatterPlotData.getMutationCount();


        double[] methylation = scatterPlotData.getDataValues(TrackType.DNA_METHYLATION.toString());
        double[] expression = scatterPlotData.getDataValues(TrackType.GENE_EXPRESSION.toString());
        double[] copyNumber = scatterPlotData.getDataValues(TrackType.COPY_NUMBER.toString());

        // check for valid data and axis selection - and if error, return null
        if (yValues == null | xValues == null | yValues.length != xValues.length)
            return;

        //classifyComboBox.setVisible(scatterPlotData.getCategories().size() > 0);

        // Note: Currently only one attribute can be selected for scatter plot,
        // therefore using the 1st series (key = index 0) to identify the attribute

        String selectedCategory = (String) classifyComboBox.getSelectedItem();
        DecimalFormat formatter = new DecimalFormat("0.00");

        XYDataModel model = new XYDataModel(selectedCategory, xAxisName, yAxisName, scatterPlotData);

        if (selectedCategory == null || ScatterPlot.isDataCategory(selectedCategory)) {
            XYSeries xySeries = new XYSeries("");
            for (int idx = 0; idx < xValues.length; ++idx) {
                // get tooltips and assign data point to series
                StringBuffer tooltip = new StringBuffer("<html>");
                tooltip.append(sampleNames[idx]);
                tooltip.append("<br>");

                for (String dn : scatterPlotData.getDataNames()) {
                    double value = scatterPlotData.getDataKeyValue(dn, idx);
                    if (!Double.isNaN(value)) {
                        tooltip.append(dn);
                        tooltip.append("=");
                        tooltip.append(formatter.format(value));
                        tooltip.append("<br>");
                    }
                }
                tooltip.append("Mutation count=");
                tooltip.append(String.valueOf(mutationCount[idx]));


                xySeries.add(idx, xValues[idx], yValues[idx], mutationCount[idx], tooltip.toString());
            }
            model.addSeries(xySeries);

        } else {
            String[] seriesNames = scatterPlotData.getAttributeCategories(selectedCategory);

            // extract sample categories for the selected attribute name
            String[] attributeValues = scatterPlotData.getSymbolValues(selectedCategory);

            // create series collection to hold xy series datasets for JFreeChart

            for (String series : seriesNames) {
                XYSeries xySeries = new XYSeries(series);
                // create plot series
                for (int idx = 0; idx < xValues.length; ++idx) {
                    // if attribute value is same as category - assign data point to the series
                    String attributeValue = attributeValues[idx];
                    if (attributeValue == null) attributeValue = "";

                    if (series.equals(attributeValue)) {

                        // get tooltips and assign data point to series
                        StringBuffer tooltip = new StringBuffer("<html>");
                        tooltip.append(sampleNames[idx]);
                        tooltip.append("<br>");
                        if (selectedCategory != null && !selectedCategory.equals("")) {
                            tooltip.append(selectedCategory);
                            tooltip.append("=");
                            tooltip.append(series);
                            tooltip.append("<br>");
                        }
                        for (String dn : scatterPlotData.getDataNames()) {
                            double value = scatterPlotData.getDataKeyValue(dn, idx);
                            if (!Double.isNaN(value)) {
                                tooltip.append(dn);
                                tooltip.append("=");
                                tooltip.append(formatter.format(value));
                                tooltip.append("<br>");
                            }
                        }
                        tooltip.append("Mutation count=");
                        tooltip.append(String.valueOf(mutationCount[idx]));

                        xySeries.add(idx, xValues[idx], yValues[idx], mutationCount[idx], tooltip.toString());
                    }
                    model.addSeries(xySeries);
                }

                // add series  dataset to series collection
            }
        }


        ScatterPlot scatterPlot = new ScatterPlot(scatterPlotData);
        scatterPlot.setModel(model);
        chartPanel.setScatterPlotModel(scatterPlot);

    }


    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        menuBar = new JMenuBar();
        fileMenu = new JMenu();
        closeMenuItem = new JMenuItem();
        commandBar = new JPanel();
        hSpacer3 = new JPanel(null);
        panel2 = new JPanel();
        label1 = new JLabel();
        xAxisComboBox = new JComboBox();
        hSpacer1 = new JPanel(null);
        panel3 = new JPanel();
        label2 = new JLabel();
        yAxisComboBox = new JComboBox();
        hSpacer2 = new JPanel(null);
        panel4 = new JPanel();
        label3 = new JLabel();
        classifyComboBox = new JComboBox();
        hSpacer4 = new JPanel(null);
        chartPanel = new ChartPanel();

        //======== this ========
        setFocusTraversalPolicyProvider(true);
        setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== menuBar ========
        {

            //======== fileMenu ========
            {
                fileMenu.setText("File");

                //---- closeMenuItem ----
                closeMenuItem.setText("Close");
                closeMenuItem.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        closeMenuItemActionPerformed(e);
                    }
                });
                fileMenu.add(closeMenuItem);
            }
            menuBar.add(fileMenu);
        }
        setJMenuBar(menuBar);

        //======== commandBar ========
        {
            commandBar.setBorder(new BevelBorder(BevelBorder.LOWERED));
            commandBar.setPreferredSize(new Dimension(420, 45));
            commandBar.setLayout(new BoxLayout(commandBar, BoxLayout.X_AXIS));
            commandBar.add(hSpacer3);

            //======== panel2 ========
            {
                panel2.setBorder(null);
                panel2.setPreferredSize(new Dimension(180, 41));
                panel2.setLayout(new BoxLayout(panel2, BoxLayout.X_AXIS));

                //---- label1 ----
                label1.setText("X: ");
                panel2.add(label1);

                //---- xAxisComboBox ----
                xAxisComboBox.setBorder(null);
                xAxisComboBox.setPreferredSize(new Dimension(120, 27));
                xAxisComboBox.setToolTipText("Set data type for X axis");
                xAxisComboBox.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        axisChanged(e);
                    }
                });
                panel2.add(xAxisComboBox);
            }
            commandBar.add(panel2);
            commandBar.add(hSpacer1);

            //======== panel3 ========
            {
                panel3.setBorder(null);
                panel3.setPreferredSize(new Dimension(180, 41));
                panel3.setLayout(new BoxLayout(panel3, BoxLayout.X_AXIS));

                //---- label2 ----
                label2.setText("Y:");
                panel3.add(label2);

                //---- yAxisComboBox ----
                yAxisComboBox.setPreferredSize(new Dimension(120, 27));
                yAxisComboBox.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        axisChanged(e);
                    }
                });
                panel3.add(yAxisComboBox);
            }
            commandBar.add(panel3);
            commandBar.add(hSpacer2);

            //======== panel4 ========
            {
                panel4.setBorder(null);
                panel4.setPreferredSize(new Dimension(250, 27));
                panel4.setLayout(new BoxLayout(panel4, BoxLayout.X_AXIS));

                //---- label3 ----
                label3.setText("Classify By:");
                panel4.add(label3);

                //---- classifyComboBox ----
                classifyComboBox.setPreferredSize(new Dimension(150, 27));
                classifyComboBox.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        attributeChanged(e);
                    }
                });
                panel4.add(classifyComboBox);
            }
            commandBar.add(panel4);
            commandBar.add(hSpacer4);
        }
        contentPane.add(commandBar, BorderLayout.NORTH);
        contentPane.add(chartPanel, BorderLayout.CENTER);
        setSize(815, 660);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JMenuBar menuBar;
    private JMenu fileMenu;
    private JMenuItem closeMenuItem;
    private JPanel commandBar;
    private JPanel hSpacer3;
    private JPanel panel2;
    private JLabel label1;
    private JComboBox xAxisComboBox;
    private JPanel hSpacer1;
    private JPanel panel3;
    private JLabel label2;
    private JComboBox yAxisComboBox;
    private JPanel hSpacer2;
    private JPanel panel4;
    private JLabel label3;
    private JComboBox classifyComboBox;
    private JPanel hSpacer4;
    private ChartPanel chartPanel;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
