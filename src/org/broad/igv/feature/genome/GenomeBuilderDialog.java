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
 * Created by JFormDesigner on Wed Apr 11 16:44:25 EDT 2012
 */

package org.broad.igv.feature.genome;

import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.igv.ui.util.MessageUtils;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

/**
 * @author Stan Diamond
 */
public class GenomeBuilderDialog extends JDialog {

    boolean canceled;
    File archiveFile;

    public GenomeBuilderDialog(Frame owner, IGV igv) {
        super(owner);
        setModal(true);
        initComponents();
        this.panel1.setIgv(igv);
    }

    public boolean isCanceled() {
        return canceled;
    }

    public File getArchiveFile() {
        return archiveFile;
    }

    private void okButtonActionPerformed(ActionEvent e) {
        boolean valid = validateFields();
        if (valid) {
            archiveFile = chooseArchiveFile();
            canceled = (archiveFile == null);
            setVisible(false);
            dispose();
        }
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        canceled = true;
        setVisible(false);
        dispose();
    }

    private File chooseArchiveFile() {
        File dir = PreferenceManager.getInstance().getLastGenomeImportDirectory();
        File initFile = new File(getGenomeId() + ".genome");
        return FileDialogUtils.chooseFile("Save .genome file", dir, initFile, FileDialog.SAVE);
    }

    private boolean validateFields() {
        StringBuffer errors = new StringBuffer();
        String id = getGenomeId();
        String name = getGenomeDisplayName();
        String fastaFile = getFastaFileName();
        if (id == null || id.length() == 0) {
            errors.append(errors.length() == 0 ? "<html>" : "<br>");
            errors.append("Unique ID is required");
        }
        if (name == null || name.length() == 0) {
            errors.append(errors.length() == 0 ? "<html>" : "<br>");
            errors.append("Descriptive name is required");
        }
        if (fastaFile == null || fastaFile.length() == 0) {
            errors.append(errors.length() == 0 ? "<html>" : "<br>");
            errors.append("FASTA file is required");
        }
        if (errors.length() == 0) {
            return true;
        } else {
            MessageUtils.showMessage(errors.toString());
            return false;
        }

    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        panel1 = new GenomeBuilderPane();
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
                contentPanel.setLayout(new BorderLayout());
                contentPanel.add(panel1, BorderLayout.CENTER);
            }
            dialogPane.add(contentPanel, BorderLayout.CENTER);

            //======== buttonBar ========
            {
                buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                buttonBar.setLayout(new GridBagLayout());
                ((GridBagLayout) buttonBar.getLayout()).columnWidths = new int[]{0, 85, 80};
                ((GridBagLayout) buttonBar.getLayout()).columnWeights = new double[]{1.0, 0.0, 0.0};

                //---- okButton ----
                okButton.setText("OK");
                okButton.addActionListener(new ActionListener() {
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
        setSize(870, 430);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private GenomeBuilderPane panel1;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables

    public String getCytobandFileName() {
        return panel1.getCytobandFileName();
    }

    public String getGeneAnnotFileName() {
        return panel1.getGeneAnnotFileName();
    }

    public String getFastaFileName() {
        return panel1.getFastaFileName();
    }

    public String getChrAliasFileName() {
        return panel1.getChrAliasFileName();
    }


    public String getGenomeDisplayName() {
        return panel1.getGenomeDisplayName();
    }

    public String getGenomeId() {
        return panel1.getGenomeId();
    }

    public String getArchiveFileName() {
        return panel1.getArchiveFileName();
    }
}
