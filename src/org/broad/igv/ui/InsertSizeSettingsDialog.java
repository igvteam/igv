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
 * Created by JFormDesigner on Thu Mar 10 19:17:05 EST 2011
 */

package org.broad.igv.ui;

import java.awt.event.*;

import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.ui.util.MessageUtils;

import java.awt.*;
import javax.swing.*;
import javax.swing.border.*;

/**
 * @author Jim Robinson
 */
public class InsertSizeSettingsDialog extends JDialog {

    private boolean isCanceled = false;

    private boolean computeIsize;
    private double minPercentile;
    private double maxPercentile;
    private int minThreshold;
    private int maxThreshold;


    public InsertSizeSettingsDialog(Frame owner, AlignmentTrack.RenderOptions options) {
        super(owner);
        initComponents();
        initValues(options);
    }

    private void initValues(AlignmentTrack.RenderOptions options) {

        computeIsize = options.isComputeIsizes();
        minPercentile = options.getMinInsertSizePercentile();
        maxPercentile = options.getMaxInsertSizePercentile();
        minThreshold = options.getMinInsertSize();
        maxThreshold = options.getMaxInsertSize();

        computeIsizeCB.setSelected(computeIsize);
        minPercentileField.setText(String.valueOf(minPercentile));
        maxPercentileField.setText(String.valueOf(maxPercentile));
        minThresholdField.setText(String.valueOf(minThreshold));
        maxThresholdField.setText(String.valueOf(maxThreshold));

        minPercentileField.setEnabled(computeIsize);
        maxPercentileField.setEnabled(computeIsize);
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        isCanceled = true;
        setVisible(false);
        dispose();
    }

    private void okButtonActionPerformed(ActionEvent e) {
        setVisible(false);
        dispose();
    }


    private void minThresholdFieldFocusLost(FocusEvent e) {
        minThresholdFieldActionPerformed(null);
    }


    private void minPercentileFieldFocusLost(FocusEvent e) {
        minPercentileFieldActionPerformed(null);
    }

    private void minThresholdFieldActionPerformed(ActionEvent e) {
        try {
            int tmp = Integer.parseInt(minThresholdField.getText());
            minThreshold = tmp;
        }
        catch (NumberFormatException ex) {
            MessageUtils.showMessage("Error: Default minimum threshold must be an integer.");
            minThresholdField.setText(String.valueOf(minThresholdField));
        }
    }


    private void maxThresholdFieldFocusLost(FocusEvent e) {
        maxThresholdFieldActionPerformed(null);
    }

    private void maxThresholdFieldActionPerformed(ActionEvent e) {
        try {
            int tmp = Integer.parseInt(maxThresholdField.getText());
            maxThreshold = tmp;
        }
        catch (NumberFormatException ex) {
            MessageUtils.showMessage("Error: Default maximum threshold must be an integer.");
            maxThresholdField.setText(String.valueOf(maxThreshold));
        }
    }


    private void maxPercentileFieldFocusLost(FocusEvent e) {
        maxPercentileFieldActionPerformed(null);
    }

    private void maxPercentileFieldActionPerformed(ActionEvent e) {
        try {
            double tmp = Double.parseDouble(maxPercentileField.getText());
            if (tmp <= 0 || tmp >= 100) throw new NumberFormatException();
            maxPercentile = tmp;
        }
        catch (NumberFormatException ex) {
            MessageUtils.showMessage("Error: Default maximum threshold must be a number between 0 and 100.");
            maxPercentileField.setText(String.valueOf(maxPercentile));
        }
    }


    private void minPercentileFieldActionPerformed(ActionEvent e) {
        try {
            double tmp = Double.parseDouble(minPercentileField.getText());
            if (tmp <= 0 || tmp >= 100) throw new NumberFormatException();
            minPercentile = tmp;
        }
        catch (NumberFormatException ex) {
            MessageUtils.showMessage("Error: Default minimum threshold must be a number between 0 and 100.");
            minPercentileField.setText(String.valueOf(minPercentile));
        }
    }

    private void computeIsizeCBActionPerformed(ActionEvent e) {
        computeIsize = computeIsizeCB.isSelected();
        minPercentileField.setEnabled(computeIsize);
        maxPercentileField.setEnabled(computeIsize);
    }


    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        label1 = new JLabel();
        panel1 = new JPanel();
        computeIsizeCB = new JCheckBox();
        label4 = new JLabel();
        label5 = new JLabel();
        minPercentileField = new JTextField();
        maxPercentileField = new JTextField();
        panel2 = new JPanel();
        label2 = new JLabel();
        label3 = new JLabel();
        minThresholdField = new JTextField();
        maxThresholdField = new JTextField();
        buttonBar = new JPanel();
        cancelButton = new JButton();
        okButton = new JButton();

        //======== this ========
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setBorder(null);
                contentPanel.setLayout(null);

                //---- label1 ----
                label1.setText("<html>These settings control color-coding of paired reads based on the inferred insert size.");
                contentPanel.add(label1);
                label1.setBounds(20, 15, 405, 60);

                //======== panel1 ========
                {
                    panel1.setBorder(new TitledBorder("Dynamic options"));
                    panel1.setLayout(null);

                    //---- computeIsizeCB ----
                    computeIsizeCB.setText("Compute thresholds");
                    computeIsizeCB.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent e) {
                            computeIsizeCBActionPerformed(e);
                        }
                    });
                    panel1.add(computeIsizeCB);
                    computeIsizeCB.setBounds(new Rectangle(new Point(10, 25), computeIsizeCB.getPreferredSize()));

                    //---- label4 ----
                    label4.setText("Minimum percentile:");
                    panel1.add(label4);
                    label4.setBounds(10, 70, 175, label4.getPreferredSize().height);

                    //---- label5 ----
                    label5.setText("Maximum percentile:");
                    panel1.add(label5);
                    label5.setBounds(10, 100, 175, 16);

                    //---- minPercentileField ----
                    minPercentileField.setEnabled(false);
                    minPercentileField.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent e) {
                            minPercentileFieldActionPerformed(e);
                        }
                    });
                    minPercentileField.addFocusListener(new FocusAdapter() {
                        @Override
                        public void focusLost(FocusEvent e) {
                            minPercentileFieldFocusLost(e);
                        }
                    });
                    panel1.add(minPercentileField);
                    minPercentileField.setBounds(220, 64, 135, minPercentileField.getPreferredSize().height);

                    //---- maxPercentileField ----
                    maxPercentileField.setEnabled(false);
                    maxPercentileField.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent e) {
                            maxPercentileFieldActionPerformed(e);
                            maxPercentileFieldActionPerformed(e);
                            maxPercentileFieldActionPerformed(e);
                        }
                    });
                    maxPercentileField.addFocusListener(new FocusAdapter() {
                        @Override
                        public void focusLost(FocusEvent e) {
                            maxPercentileFieldFocusLost(e);
                        }
                    });
                    panel1.add(maxPercentileField);
                    maxPercentileField.setBounds(220, 94, 135, 28);

                    { // compute preferred size
                        Dimension preferredSize = new Dimension();
                        for(int i = 0; i < panel1.getComponentCount(); i++) {
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
                contentPanel.add(panel1);
                panel1.setBounds(25, 210, 430, 160);

                //======== panel2 ========
                {
                    panel2.setBorder(new TitledBorder("Defaults"));
                    panel2.setLayout(null);

                    //---- label2 ----
                    label2.setText("Default minimum threshold: ");
                    panel2.add(label2);
                    label2.setBounds(10, 30, 205, label2.getPreferredSize().height);

                    //---- label3 ----
                    label3.setText("Default maximum threshold: ");
                    panel2.add(label3);
                    label3.setBounds(10, 60, 185, 16);

                    //---- minThresholdField ----
                    minThresholdField.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent e) {
                            minThresholdFieldActionPerformed(e);
                        }
                    });
                    minThresholdField.addFocusListener(new FocusAdapter() {
                        @Override
                        public void focusLost(FocusEvent e) {
                            minThresholdFieldFocusLost(e);
                        }
                    });
                    panel2.add(minThresholdField);
                    minThresholdField.setBounds(215, 24, 145, minThresholdField.getPreferredSize().height);

                    //---- maxThresholdField ----
                    maxThresholdField.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent e) {
                            maxThresholdFieldActionPerformed(e);
                        }
                    });
                    maxThresholdField.addFocusListener(new FocusAdapter() {
                        @Override
                        public void focusLost(FocusEvent e) {
                            maxThresholdFieldFocusLost(e);
                        }
                    });
                    panel2.add(maxThresholdField);
                    maxThresholdField.setBounds(215, 54, 145, 28);

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
                contentPanel.add(panel2);
                panel2.setBounds(25, 85, 430, 110);

                { // compute preferred size
                    Dimension preferredSize = new Dimension();
                    for(int i = 0; i < contentPanel.getComponentCount(); i++) {
                        Rectangle bounds = contentPanel.getComponent(i).getBounds();
                        preferredSize.width = Math.max(bounds.x + bounds.width, preferredSize.width);
                        preferredSize.height = Math.max(bounds.y + bounds.height, preferredSize.height);
                    }
                    Insets insets = contentPanel.getInsets();
                    preferredSize.width += insets.right;
                    preferredSize.height += insets.bottom;
                    contentPanel.setMinimumSize(preferredSize);
                    contentPanel.setPreferredSize(preferredSize);
                }
            }
            dialogPane.add(contentPanel, BorderLayout.CENTER);

            //======== buttonBar ========
            {
                buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                buttonBar.setLayout(new FlowLayout(FlowLayout.RIGHT));

                //---- cancelButton ----
                cancelButton.setText("Cancel");
                cancelButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        cancelButtonActionPerformed(e);
                        cancelButtonActionPerformed(e);
                    }
                });
                buttonBar.add(cancelButton);

                //---- okButton ----
                okButton.setText("OK");
                okButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        okButtonActionPerformed(e);
                    }
                });
                buttonBar.add(okButton);
            }
            dialogPane.add(buttonBar, BorderLayout.SOUTH);
        }
        contentPane.add(dialogPane, BorderLayout.CENTER);
        pack();
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JLabel label1;
    private JPanel panel1;
    private JCheckBox computeIsizeCB;
    private JLabel label4;
    private JLabel label5;
    private JTextField minPercentileField;
    private JTextField maxPercentileField;
    private JPanel panel2;
    private JLabel label2;
    private JLabel label3;
    private JTextField minThresholdField;
    private JTextField maxThresholdField;
    private JPanel buttonBar;
    private JButton cancelButton;
    private JButton okButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables

    public boolean isCanceled() {
        return isCanceled;
    }

    public boolean isComputeIsize() {
        return computeIsize;
    }

    public double getMinPercentile() {
        return minPercentile;
    }

    public double getMaxPercentile() {
        return maxPercentile;
    }

    public int getMinThreshold() {
        return minThreshold;
    }

    public int getMaxThreshold() {
        return maxThreshold;
    }
}
