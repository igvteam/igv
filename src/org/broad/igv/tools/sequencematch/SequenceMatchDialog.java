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
 * Created by JFormDesigner on Tue Jan 29 15:22:45 EST 2013
 */

package org.broad.igv.tools.sequencematch;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * Dialog so the user can enter a pattern that can be used
 * to search a nucleotide sequence.
 *
 * @author Jacob Silterra
 */
public class SequenceMatchDialog extends JDialog {

    private String trackName;
    private String inputPattern;

    public SequenceMatchDialog(Frame owner) {
        super(owner);
        initComponents();
    }

    public String getInputPattern() {
        return inputPattern;
    }

    public String getTrackName() {
        return trackName;
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        this.setVisible(false);
    }

    private void okButtonActionPerformed(ActionEvent e) {

        this.trackName = null;
        this.inputPattern = null;

        String strPattern = patternField.getText().toUpperCase();

        boolean patternIsValid = checkPatternValid(strPattern);
        if(!patternIsValid){
            // TODO warn user
        }else{
            this.trackName = nameField.getText();
            this.inputPattern = strPattern;
            this.setVisible(false);
        }
    }

    private boolean checkPatternValid(String strPattern) {
        for(int ii=0; ii < strPattern.length(); ii++){
            String c = strPattern.substring(ii, ii + 1);
            if(!SequenceMatchSource.isValidString(c)){
                return false;
            }
        }
        return true;
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        textArea1 = new JTextArea();
        label1 = new JLabel();
        nameField = new JTextField();
        label2 = new JLabel();
        label4 = new JLabel();
        patternField = new JTextField();
        buttonBar = new JPanel();
        okButton = new JButton();
        cancelButton = new JButton();

        //======== this ========
        setModal(true);
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setAlignmentX(0.0F);
                contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.Y_AXIS));

                //---- textArea1 ----
                textArea1.setText("Enter the pattern to search for, using ACGT and/or IUPAC ambiguity codes. Matches in the current genome sequence will be displayed in a new track.\n");
                textArea1.setEditable(false);
                textArea1.setBackground(new Color(238, 238, 238));
                textArea1.setLineWrap(true);
                contentPanel.add(textArea1);

                //---- label1 ----
                label1.setText("Track Name:");
                label1.setLabelFor(nameField);
                label1.setHorizontalTextPosition(SwingConstants.LEFT);
                label1.setHorizontalAlignment(SwingConstants.LEFT);
                label1.setMaximumSize(new Dimension(374, 16));
                label1.setPreferredSize(new Dimension(374, 16));
                label1.setAlignmentX(1.0F);
                contentPanel.add(label1);
                contentPanel.add(nameField);

                //---- label2 ----
                label2.setText("Pattern:");
                label2.setLabelFor(patternField);
                label2.setHorizontalTextPosition(SwingConstants.LEFT);
                label2.setHorizontalAlignment(SwingConstants.LEFT);
                label2.setAlignmentX(1.0F);
                label2.setMaximumSize(new Dimension(374, 16));
                label2.setPreferredSize(new Dimension(374, 16));
                contentPanel.add(label2);
                contentPanel.add(label4);
                contentPanel.add(patternField);
            }
            dialogPane.add(contentPanel, BorderLayout.NORTH);

            //======== buttonBar ========
            {
                buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                buttonBar.setLayout(new GridBagLayout());
                ((GridBagLayout)buttonBar.getLayout()).columnWidths = new int[] {0, 85, 80};
                ((GridBagLayout)buttonBar.getLayout()).columnWeights = new double[] {1.0, 0.0, 0.0};

                //---- okButton ----
                okButton.setText("OK");
                okButton.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        okButtonActionPerformed(e);
                    }
                });
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
            dialogPane.add(buttonBar, BorderLayout.SOUTH);
        }
        contentPane.add(dialogPane, BorderLayout.CENTER);
        setSize(400, 300);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JTextArea textArea1;
    private JLabel label1;
    private JTextField nameField;
    private JLabel label2;
    private JLabel label4;
    private JTextField patternField;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
