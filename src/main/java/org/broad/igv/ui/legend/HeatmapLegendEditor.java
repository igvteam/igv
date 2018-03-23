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
 * Created by JFormDesigner on Thu Jun 16 11:12:56 EDT 2011
 */

package org.broad.igv.ui.legend;

import java.awt.*;
import java.awt.Component;
import java.awt.event.*;
import javax.swing.*;

import org.broad.igv.renderer.ColorScale;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.color.ColorChooserPanel;
import org.jdesktop.layout.GroupLayout;
import org.jdesktop.layout.LayoutStyle;

/**
 * @author Stan Diamond
 */
public class HeatmapLegendEditor extends JDialog {
    private boolean canceled = true;
    private ContinuousColorScale colorScheme;
    private TrackType type;

    /**
     * Creates new form HeatmapLegendEditor2
     */
    public HeatmapLegendEditor(java.awt.Frame parent, boolean modal, TrackType type, ColorScale colorScheme) {
        super(parent, modal);
        this.colorScheme = (ContinuousColorScale) colorScheme;
        this.type = type;
        initComponents();
        initValues();
        this.setLocationRelativeTo(parent);
        this.getRootPane().setDefaultButton(okButton);
    }

    private void initValues() {
        doubleGradientCheckbox.setSelected(colorScheme.isUseDoubleGradient());
        negRangeStart.setText(String.valueOf(colorScheme.getNegStart()));
        negRangeEnd.setText(String.valueOf(getColorScheme().getMinimum()));
        posRangeStart.setText(String.valueOf(colorScheme.getPosStart()));
        posRangeEnd.setText(String.valueOf(colorScheme.getMaximum()));
        minColor.setSelectedColor(colorScheme.getMinColor());
        maxColor.setSelectedColor(colorScheme.getMaxColor());

        // Single gradient color schems might have a null mid color.  Default
        // to white in that case, a non-null color is required.
        Color mc = colorScheme.getMidColor();
        midColor.setSelectedColor(mc == null ? Color.white : mc);

        initDoubleGradientState();

    }

    private void initDoubleGradientState() {
        final boolean doubleGradient = doubleGradientCheckbox.isSelected();
        negRangePanel.setVisible(doubleGradient);
        midColorLabel.setVisible(doubleGradient);
        midColor.setVisible(doubleGradient);
        posRangeLabel.setText(doubleGradient ? "Positive Range " : "Range");
    }

    private boolean updateValues() {

        try {
            double negStart = 0;
            double negEnd = 0;
            double posStart = Double.parseDouble(posRangeStart.getText());
            double posEnd = Double.parseDouble(posRangeEnd.getText());
            negStart = Double.parseDouble(negRangeStart.getText());
            negEnd = Double.parseDouble(negRangeEnd.getText());


            colorScheme = new ContinuousColorScale(
                    Math.max(negStart, negEnd),
                    Math.min(negStart, negEnd),
                    Math.min(posStart, posEnd),
                    Math.max(posStart, posEnd),
                    minColor.getSelectedColor(),
                    midColor.getSelectedColor(),
                    maxColor.getSelectedColor());

            return true;

        } catch (NumberFormatException numberFormatException) {
            JOptionPane.showMessageDialog(this, "Limit fields must be numeric.", "Error", JOptionPane.ERROR_MESSAGE);
            return false;
        }
    }


    public ContinuousColorScale getColorScheme() {
        return colorScheme;
    }


    private void okButtonActionPerformed(ActionEvent e) {
        canceled = false;
        if (updateValues()) {
            setVisible(false);
        }

    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        canceled = true;
        setVisible(false);
    }

    private void negRangeEndActionPerformed(ActionEvent e) {
        // TODO add your code here
    }

    private void doubleGradientCheckboxActionPerformed(ActionEvent e) {
         initDoubleGradientState();
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        jPanel1 = new JPanel();
        midColorLabel = new JLabel();
        jLabel3 = new JLabel();
        minColorLabel = new JLabel();
        minColor = new ColorChooserPanel();
        midColor = new ColorChooserPanel();
        maxColor = new ColorChooserPanel();
        okButton = new JButton();
        cancelButton = new JButton();
        negRangePanel = new JPanel();
        negRangeLabel = new JLabel();
        negRangeStart = new JTextField();
        negRangeToLabel = new JLabel();
        negRangeEnd = new JTextField();
        doubleGradientCheckbox = new JCheckBox();
        posRangePanel = new JPanel();
        posRangeLabel = new JLabel();
        posRangeStart = new JTextField();
        posRangeToLabel = new JLabel();
        posRangeEnd = new JTextField();

        //======== this ========
        setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        setResizable(false);
        Container contentPane = getContentPane();

        //======== jPanel1 ========
        {

            //---- midColorLabel ----
            midColorLabel.setText("Midpoint Color");

            //---- jLabel3 ----
            jLabel3.setText("Maximum Color");

            //---- minColorLabel ----
            minColorLabel.setText("Minimum Color");


            GroupLayout jPanel1Layout = new GroupLayout(jPanel1);
            jPanel1.setLayout(jPanel1Layout);
            jPanel1Layout.setHorizontalGroup(
                    jPanel1Layout.createParallelGroup()
                            .add(jPanel1Layout.createSequentialGroup()
                            .add(jPanel1Layout.createParallelGroup()
                                    .add(jPanel1Layout.createSequentialGroup()
                                            .add(jLabel3)
                                            .addPreferredGap(LayoutStyle.RELATED, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                            .add(maxColor, GroupLayout.PREFERRED_SIZE, 50, GroupLayout.PREFERRED_SIZE))
                                    .add(jPanel1Layout.createSequentialGroup()
                                            .add(minColorLabel)
                                            .addPreferredGap(LayoutStyle.RELATED, 9, Short.MAX_VALUE)
                                            .add(minColor, GroupLayout.PREFERRED_SIZE, 75, GroupLayout.PREFERRED_SIZE))
                                    .add(GroupLayout.TRAILING, jPanel1Layout.createSequentialGroup()
                                    .add(midColorLabel)
                                    .addPreferredGap(LayoutStyle.RELATED, 12, Short.MAX_VALUE)
                                    .add(midColor, GroupLayout.PREFERRED_SIZE, 50, GroupLayout.PREFERRED_SIZE)))
                            .addContainerGap())
            );
            jPanel1Layout.linkSize(new Component[]{maxColor, midColor, minColor}, GroupLayout.HORIZONTAL);
            jPanel1Layout.setVerticalGroup(
                    jPanel1Layout.createParallelGroup()
                            .add(jPanel1Layout.createSequentialGroup()
                            .add(jPanel1Layout.createParallelGroup()
                                    .add(jPanel1Layout.createSequentialGroup()
                                            .add(minColorLabel)
                                            .addPreferredGap(LayoutStyle.UNRELATED)
                                            .add(midColorLabel))
                                    .add(jPanel1Layout.createSequentialGroup()
                                    .add(minColor, GroupLayout.PREFERRED_SIZE, 25, GroupLayout.PREFERRED_SIZE)
                                    .add(14, 14, 14)
                                    .add(midColor, GroupLayout.PREFERRED_SIZE, 25, GroupLayout.PREFERRED_SIZE)))
                            .add(12, 12, 12)
                            .add(jPanel1Layout.createParallelGroup()
                            .add(GroupLayout.TRAILING, jPanel1Layout.createSequentialGroup()
                                    .add(jLabel3)
                                    .add(20, 20, 20))
                            .add(jPanel1Layout.createSequentialGroup()
                            .add(maxColor, GroupLayout.PREFERRED_SIZE, 25, GroupLayout.PREFERRED_SIZE)
                            .addContainerGap())))
            );
            jPanel1Layout.linkSize(new Component[]{jLabel3, maxColor}, GroupLayout.VERTICAL);
            jPanel1Layout.linkSize(new Component[]{midColor, midColorLabel}, GroupLayout.VERTICAL);
            jPanel1Layout.linkSize(new Component[]{minColor, minColorLabel}, GroupLayout.VERTICAL);
        }

        //---- okButton ----
        okButton.setText("OK");
        okButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                okButtonActionPerformed(e);
            }
        });

        //---- cancelButton ----
        cancelButton.setText("Cancel");
        cancelButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                cancelButtonActionPerformed(e);
            }
        });

        //======== negRangePanel ========
        {

            //---- negRangeLabel ----
            negRangeLabel.setText("Negative Range: ");

            //---- negRangeStart ----
            negRangeStart.setText("-0.1");

            //---- negRangeToLabel ----
            negRangeToLabel.setText("To:");

            //---- negRangeEnd ----
            negRangeEnd.setText("-1.5");
            negRangeEnd.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    negRangeEndActionPerformed(e);
                }
            });

            GroupLayout negRangePanelLayout = new GroupLayout(negRangePanel);
            negRangePanel.setLayout(negRangePanelLayout);
            negRangePanelLayout.setHorizontalGroup(
                    negRangePanelLayout.createParallelGroup()
                            .add(negRangePanelLayout.createSequentialGroup()
                            .add(negRangeLabel, GroupLayout.PREFERRED_SIZE, 115, GroupLayout.PREFERRED_SIZE)
                            .addPreferredGap(LayoutStyle.RELATED)
                            .add(negRangeStart, GroupLayout.PREFERRED_SIZE, 59, GroupLayout.PREFERRED_SIZE)
                            .addPreferredGap(LayoutStyle.UNRELATED)
                            .add(negRangeToLabel, GroupLayout.PREFERRED_SIZE, 30, GroupLayout.PREFERRED_SIZE)
                            .addPreferredGap(LayoutStyle.RELATED)
                            .add(negRangeEnd, GroupLayout.PREFERRED_SIZE, 59, GroupLayout.PREFERRED_SIZE))
            );
            negRangePanelLayout.setVerticalGroup(
                    negRangePanelLayout.createParallelGroup()
                            .add(negRangePanelLayout.createSequentialGroup()
                            .add(negRangePanelLayout.createParallelGroup()
                                    .add(negRangeLabel, GroupLayout.PREFERRED_SIZE, 22, GroupLayout.PREFERRED_SIZE)
                                    .add(negRangePanelLayout.createParallelGroup(GroupLayout.BASELINE)
                                    .add(negRangeStart, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                    .add(negRangeToLabel, GroupLayout.PREFERRED_SIZE, 22, GroupLayout.PREFERRED_SIZE)
                                    .add(negRangeEnd, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)))
                            .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
            );
        }

        //---- doubleGradientCheckbox ----
        doubleGradientCheckbox.setText("Use Double Gradient");
        doubleGradientCheckbox.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                doubleGradientCheckboxActionPerformed(e);
            }
        });

        //======== posRangePanel ========
        {

            //---- posRangeLabel ----
            posRangeLabel.setText("Positive Range: ");

            //---- posRangeStart ----
            posRangeStart.setText("-0.1");

            //---- posRangeToLabel ----
            posRangeToLabel.setText("To:");

            //---- posRangeEnd ----
            posRangeEnd.setText("-1.5");
            posRangeEnd.setMaximumSize(new Dimension(36, 22));
            posRangeEnd.setMinimumSize(new Dimension(36, 22));

            GroupLayout posRangePanelLayout = new GroupLayout(posRangePanel);
            posRangePanel.setLayout(posRangePanelLayout);
            posRangePanelLayout.setHorizontalGroup(
                    posRangePanelLayout.createParallelGroup()
                            .add(posRangePanelLayout.createSequentialGroup()
                            .add(posRangeLabel, GroupLayout.PREFERRED_SIZE, 115, GroupLayout.PREFERRED_SIZE)
                            .addPreferredGap(LayoutStyle.RELATED)
                            .add(posRangeStart, GroupLayout.PREFERRED_SIZE, 59, GroupLayout.PREFERRED_SIZE)
                            .addPreferredGap(LayoutStyle.UNRELATED)
                            .add(posRangeToLabel, GroupLayout.PREFERRED_SIZE, 30, GroupLayout.PREFERRED_SIZE)
                            .addPreferredGap(LayoutStyle.RELATED)
                            .add(posRangeEnd, GroupLayout.PREFERRED_SIZE, 59, GroupLayout.PREFERRED_SIZE))
            );
            posRangePanelLayout.setVerticalGroup(
                    posRangePanelLayout.createParallelGroup()
                            .add(posRangePanelLayout.createSequentialGroup()
                            .add(posRangePanelLayout.createParallelGroup()
                                    .add(posRangeLabel, GroupLayout.PREFERRED_SIZE, 22, GroupLayout.PREFERRED_SIZE)
                                    .add(posRangePanelLayout.createParallelGroup(GroupLayout.BASELINE)
                                    .add(posRangeStart, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                    .add(posRangeToLabel, GroupLayout.PREFERRED_SIZE, 22, GroupLayout.PREFERRED_SIZE)
                                    .add(posRangeEnd, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)))
                            .addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
            );
        }

        GroupLayout contentPaneLayout = new GroupLayout(contentPane);
        contentPane.setLayout(contentPaneLayout);
        contentPaneLayout.setHorizontalGroup(
                contentPaneLayout.createParallelGroup()
                        .add(contentPaneLayout.createSequentialGroup()
                                .add(35, 35, 35)
                                .add(contentPaneLayout.createParallelGroup(GroupLayout.LEADING, false)
                                        .add(doubleGradientCheckbox, GroupLayout.PREFERRED_SIZE, 175, GroupLayout.PREFERRED_SIZE)
                                        .add(jPanel1, GroupLayout.PREFERRED_SIZE, 184, GroupLayout.PREFERRED_SIZE)
                                        .add(negRangePanel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                        .add(GroupLayout.TRAILING, posRangePanel, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                                .addContainerGap(101, Short.MAX_VALUE))
                        .add(GroupLayout.TRAILING, contentPaneLayout.createSequentialGroup()
                        .addContainerGap(124, Short.MAX_VALUE)
                        .add(okButton)
                        .addPreferredGap(LayoutStyle.RELATED)
                        .add(cancelButton)
                        .add(132, 132, 132))
        );
        contentPaneLayout.setVerticalGroup(
                contentPaneLayout.createParallelGroup()
                        .add(contentPaneLayout.createSequentialGroup()
                        .add(52, 52, 52)
                        .add(doubleGradientCheckbox)
                        .add(18, 18, 18)
                        .add(jPanel1, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(LayoutStyle.RELATED)
                        .add(negRangePanel, GroupLayout.PREFERRED_SIZE, 30, GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(LayoutStyle.UNRELATED)
                        .add(posRangePanel, GroupLayout.PREFERRED_SIZE, 30, GroupLayout.PREFERRED_SIZE)
                        .add(18, 18, 18)
                        .add(contentPaneLayout.createParallelGroup(GroupLayout.BASELINE)
                                .add(okButton)
                                .add(cancelButton))
                        .add(34, 34, 34))
        );
        setSize(425, 405);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel jPanel1;
    private JLabel midColorLabel;
    private JLabel jLabel3;
    private JLabel minColorLabel;
    private ColorChooserPanel minColor;
    private ColorChooserPanel midColor;
    private ColorChooserPanel maxColor;
    private JButton okButton;
    private JButton cancelButton;
    private JPanel negRangePanel;
    private JLabel negRangeLabel;
    private JTextField negRangeStart;
    private JLabel negRangeToLabel;
    private JTextField negRangeEnd;
    private JCheckBox doubleGradientCheckbox;
    private JPanel posRangePanel;
    private JLabel posRangeLabel;
    private JTextField posRangeStart;
    private JLabel posRangeToLabel;
    private JTextField posRangeEnd;
    // JFormDesigner - End of variables declaration  //GEN-END:variables

    public boolean isCanceled() {
        return canceled;
    }
}
