/*
 * Created by JFormDesigner on Fri Oct 23 16:50:52 PDT 2015
 */

package org.broad.igv.ui.util;

import org.broad.igv.Globals;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * @author James Robinson
 */
public class LoadFromURLDialog extends org.broad.igv.ui.IGVDialog {

    boolean canceled = false;
    String fileURL;
    String indexURL;
    JTextField fileField;
    JTextField indexField;

    public LoadFromURLDialog(Frame owner, boolean isHtsget) {
        super(owner, isHtsget ? "htsget URL" : "Load from URL");
        initComponents(isHtsget);
    }

    private void okButtonActionPerformed(ActionEvent e) {
        canceled = false;

        fileURL = fileField.getText().trim();
        indexURL = indexField.getText().trim();
        if (indexURL.length() == 0) {
            indexURL = null;
        }
        setVisible(false);
        dispose();
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        canceled = true;
        setVisible(false);
        dispose();
    }


    public boolean isCanceled() {
        return canceled;
    }

    public List<String> getFileURLs() {
        return splitOnWhiteSpace(fileURL);
    }

    private static List<String> splitOnWhiteSpace(String string) {
        if (string != null && !string.isBlank()) {
            String[] inputs = Globals.whitespacePattern.split(string.trim());
            return Arrays.asList(inputs);
        }
        return Collections.emptyList();
    }

    public List<String> getIndexURLs() {
        return splitOnWhiteSpace(indexURL);
    }

    private void initComponents(boolean isHtsget) {
        // Panels
        JPanel dialogPane = new JPanel(new BorderLayout(0, 12));
        dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
        JPanel contentPanel = new JPanel(new GridBagLayout());
        JPanel buttonBar = new JPanel(new GridBagLayout());
        buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));

        // Components
        fileField = new JTextField(50); // Set preferred width
        indexField = new JTextField(50); // Set preferred width
        JButton okButton = new JButton("OK");
        JButton cancelButton = new JButton("Cancel");

        //======== Content Panel ========
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(0, 0, 5, 5);

        // Row 0: File URL
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.weightx = 0.0;
        contentPanel.add(new JLabel("File URL:"), gbc);

        gbc.gridx = 1;
        gbc.weightx = 1.0;
        gbc.insets = new Insets(0, 0, 5, 0);
        contentPanel.add(fileField, gbc);

        if (!isHtsget) {
            // Row 1: Help Text
            gbc.gridy = 1;
            gbc.gridx = 0;
            gbc.gridwidth = 2;
            gbc.weightx = 0.0;
            gbc.insets = new Insets(0, 0, 5, 0);
            contentPanel.add(new JLabel("<html><i>Specify url to an index file. <b>Required for BAM and indexed files</b></i>"), gbc);

            // Row 2: Index URL
            gbc.gridy = 2;
            gbc.gridwidth = 1; // Reset gridwidth
            gbc.gridx = 0;
            gbc.weightx = 0.0;
            gbc.insets = new Insets(0, 0, 0, 5);
            contentPanel.add(new JLabel("Index URL:"), gbc);

            gbc.gridx = 1;
            gbc.weightx = 1.0;
            gbc.insets = new Insets(0, 0, 0, 0);
            contentPanel.add(indexField, gbc);
        }

        //======== Button Bar ========
        GridBagConstraints gbcButtons = new GridBagConstraints();
        gbcButtons.gridx = 0;
        gbcButtons.weightx = 1.0; // Pushes buttons to the right
        buttonBar.add(Box.createHorizontalGlue(), gbcButtons);

        gbcButtons.gridx = 1;
        gbcButtons.weightx = 0.0;
        buttonBar.add(okButton, gbcButtons);

        gbcButtons.gridx = 2;
        buttonBar.add(cancelButton, gbcButtons);

        //======== Assemble Dialog ========
        dialogPane.add(contentPanel, BorderLayout.CENTER);
        dialogPane.add(buttonBar, BorderLayout.SOUTH);
        setContentPane(dialogPane);

        // Action Listeners
        okButton.addActionListener(this::okButtonActionPerformed);
        cancelButton.addActionListener(this::cancelButtonActionPerformed);
        fileField.addActionListener(this::okButtonActionPerformed);
        indexField.addActionListener(this::okButtonActionPerformed);

        // Final setup
        setModalityType(Dialog.ModalityType.APPLICATION_MODAL);
        pack();
        setLocationRelativeTo(getOwner());
        getRootPane().setDefaultButton(okButton);
    }

}
