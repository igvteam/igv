/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

/*
 * Created by JFormDesigner on Mon Nov 26 16:20:15 EST 2012
 */

package org.broad.igv.variant;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;

import org.broadinstitute.sting.gatk.walkers.na12878kb.TruthStatus;

/**
 * @author User #2
 */
public class VariantReviewDialog extends JDialog {
    public VariantReviewDialog(Frame owner, String sample, Variant variant) {
        super(owner);
        initComponents();

        truthField.setModel(new DefaultComboBoxModel(TruthStatus.values()));
        initComponentData(sample, variant);
    }

    private void initComponentData(String sample, Variant variant) {
        String uname = System.getProperty("user.name", "unknown");
        callsetField.setText(uname);

        chrField.setText(variant.getChr());
        startField.setText("" + variant.getStart());
        stopField.setText("" + variant.getEnd());

        Genotype genotype = variant.getGenotype(sample);
        String genoString = genotype.getGenotypeString();
        genotypeField.setText(genoString);

        String type = variant.getType();
        System.out.println("type: " + type);

        truthField.getModel().setSelectedItem(TruthStatus.UNKNOWN);

        validate();
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        setVisible(false);
    }


    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        panel1 = new JPanel();
        label4 = new JLabel();
        callsetField = new JTextField();
        panel2 = new JPanel();
        label5 = new JLabel();
        truthField = new JComboBox();
        panel3 = new JPanel();
        label6 = new JLabel();
        genotypeField = new JLabel();
        panel4 = new JPanel();
        label7 = new JLabel();
        chrField = new JLabel();
        panel5 = new JPanel();
        label8 = new JLabel();
        startField = new JLabel();
        panel6 = new JPanel();
        label9 = new JLabel();
        stopField = new JLabel();
        buttonBar = new JPanel();
        okButton = new JButton();
        cancelButton = new JButton();

        //======== this ========
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.X_AXIS));

                //======== panel1 ========
                {
                    panel1.setLayout(new BoxLayout(panel1, BoxLayout.Y_AXIS));

                    //---- label4 ----
                    label4.setText("Callset:");
                    label4.setHorizontalAlignment(SwingConstants.LEFT);
                    label4.setMaximumSize(new Dimension(80, 16));
                    panel1.add(label4);
                    panel1.add(callsetField);
                }
                contentPanel.add(panel1);

                //======== panel2 ========
                {
                    panel2.setLayout(new BoxLayout(panel2, BoxLayout.Y_AXIS));

                    //---- label5 ----
                    label5.setText("Truth");
                    label5.setHorizontalAlignment(SwingConstants.LEFT);
                    label5.setMaximumSize(new Dimension(80, 16));
                    label5.setLabelFor(truthField);
                    panel2.add(label5);
                    panel2.add(truthField);
                }
                contentPanel.add(panel2);

                //======== panel3 ========
                {
                    panel3.setLayout(new BoxLayout(panel3, BoxLayout.Y_AXIS));

                    //---- label6 ----
                    label6.setText("Genotype");
                    label6.setHorizontalAlignment(SwingConstants.LEFT);
                    label6.setMaximumSize(new Dimension(80, 16));
                    panel3.add(label6);

                    //---- genotypeField ----
                    genotypeField.setText("text");
                    panel3.add(genotypeField);
                }
                contentPanel.add(panel3);

                //======== panel4 ========
                {
                    panel4.setLayout(new BoxLayout(panel4, BoxLayout.Y_AXIS));

                    //---- label7 ----
                    label7.setText("Chr");
                    label7.setHorizontalAlignment(SwingConstants.LEFT);
                    label7.setMaximumSize(new Dimension(80, 16));
                    panel4.add(label7);
                    panel4.add(chrField);
                }
                contentPanel.add(panel4);

                //======== panel5 ========
                {
                    panel5.setLayout(new BoxLayout(panel5, BoxLayout.Y_AXIS));

                    //---- label8 ----
                    label8.setText("Start");
                    label8.setHorizontalAlignment(SwingConstants.LEFT);
                    label8.setMaximumSize(new Dimension(80, 16));
                    panel5.add(label8);
                    panel5.add(startField);
                }
                contentPanel.add(panel5);

                //======== panel6 ========
                {
                    panel6.setLayout(new BoxLayout(panel6, BoxLayout.Y_AXIS));

                    //---- label9 ----
                    label9.setText("Stop");
                    label9.setHorizontalAlignment(SwingConstants.LEFT);
                    label9.setMaximumSize(new Dimension(80, 16));
                    panel6.add(label9);
                    panel6.add(stopField);
                }
                contentPanel.add(panel6);
            }
            dialogPane.add(contentPanel, BorderLayout.CENTER);

            //======== buttonBar ========
            {
                buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                buttonBar.setLayout(new GridBagLayout());
                ((GridBagLayout)buttonBar.getLayout()).columnWidths = new int[] {0, 85, 80};
                ((GridBagLayout)buttonBar.getLayout()).columnWeights = new double[] {1.0, 0.0, 0.0};

                //---- okButton ----
                okButton.setText("OK");
                buttonBar.add(okButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 0, 5), 0, 0));

                //---- cancelButton ----
                cancelButton.setText("Cancel");
                cancelButton.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        cancelButtonActionPerformed(e);
                    }
                });
                buttonBar.add(cancelButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 0, 0), 0, 0));
            }
            dialogPane.add(buttonBar, BorderLayout.PAGE_END);
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
    private JPanel panel1;
    private JLabel label4;
    private JTextField callsetField;
    private JPanel panel2;
    private JLabel label5;
    private JComboBox truthField;
    private JPanel panel3;
    private JLabel label6;
    private JLabel genotypeField;
    private JPanel panel4;
    private JLabel label7;
    private JLabel chrField;
    private JPanel panel5;
    private JLabel label8;
    private JLabel startField;
    private JPanel panel6;
    private JLabel label9;
    private JLabel stopField;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
