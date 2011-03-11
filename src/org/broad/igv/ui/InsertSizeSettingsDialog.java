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
 * Created by JFormDesigner on Thu Mar 10 19:17:05 EST 2011
 */

package org.broad.igv.ui;

import java.awt.*;
import javax.swing.*;
import javax.swing.border.*;

/**
 * @author Jim Robinson
 */
public class InsertSizeSettingsDialog extends JDialog {

    public InsertSizeSettingsDialog(Frame owner) {
        super(owner);
        initComponents();
    }

    public InsertSizeSettingsDialog(Dialog owner) {
        super(owner);
        initComponents();
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        label1 = new JLabel();
        panel1 = new JPanel();
        checkBox1 = new JCheckBox();
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

                    //---- checkBox1 ----
                    checkBox1.setText("Compute thresholds dynamically");
                    panel1.add(checkBox1);
                    checkBox1.setBounds(new Rectangle(new Point(10, 25), checkBox1.getPreferredSize()));

                    //---- label4 ----
                    label4.setText("Minimum percentile:");
                    panel1.add(label4);
                    label4.setBounds(10, 70, 175, label4.getPreferredSize().height);

                    //---- label5 ----
                    label5.setText("Maximum percentile:");
                    panel1.add(label5);
                    label5.setBounds(10, 100, 175, 16);
                    panel1.add(minPercentileField);
                    minPercentileField.setBounds(220, 64, 135, minPercentileField.getPreferredSize().height);
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
                    panel2.add(minThresholdField);
                    minThresholdField.setBounds(215, 24, 145, minThresholdField.getPreferredSize().height);
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
                buttonBar.add(cancelButton);

                //---- okButton ----
                okButton.setText("OK");
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
    private JCheckBox checkBox1;
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
}
