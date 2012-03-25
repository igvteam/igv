/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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
package org.broad.igv.feature.genome;

import java.awt.*;
import java.awt.event.*;
import javax.swing.border.*;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.FileDialogUtils;

import javax.swing.*;
import java.io.File;
import java.util.Collection;

/**
 * @author eflakes
 */
public class GenomeBuilderPane extends javax.swing.JPanel {

    private static Logger logger = Logger.getLogger(GenomeBuilderPane.class);
    private String genomeArchiveLocation;
    private String genomeFilename;
    GenomeImporter importer;
    IGV igv;

    /**
     * Creates new form GenomeBuilder
     */
    public GenomeBuilderPane(IGV igv) {
        this.igv = igv;
        initComponents();
        importer = new GenomeImporter();
    }

    public String getCytobandFileName() {
        String cytobandFile = cytobandFileTextField.getText();
        if (cytobandFile != null && cytobandFile.trim().equals("")) {
            cytobandFile = null;
        }
        return cytobandFile;
    }

    public String getRefFlatFileName() {
        String refFlatFile = refFlatFileTextField.getText();
        if (refFlatFile != null && refFlatFile.trim().equals("")) {
            refFlatFile = null;
        }
        return refFlatFile;
    }

    public String getFastaFileName() {
        String fastaFile = fastaFileTextField.getText();
        if (fastaFile != null && fastaFile.trim().equals("")) {
            fastaFile = null;
        }
        return fastaFile;
    }

    public String getChrAliasFileName() {
        String chrAliasFile = chrAliasField.getText();
        if (chrAliasFile != null && chrAliasFile.trim().equals("")) {
            chrAliasFile = null;
        }
        return chrAliasFile;
    }


    public String getGenomeId() {
        return idField.getText();
    }

    public String getSequenceURL() {
        return sequenceURLField.getText();
    }


    public String getGenomeDisplayName() {
        String name = genomeDisplayNameTextField.getText();
        if (name != null && name.trim().equals("")) {
            name = null;
        } else {
            name = name.trim();
        }
        return name;
    }

    public String getGenomeArchiveLocation() {
        if (genomeArchiveLocation != null && genomeArchiveLocation.trim().equals("")) {
            genomeArchiveLocation = null;
        }
        return genomeArchiveLocation;
    }

    public String getArchiveFileName() {

        if (genomeFilename == null) {
            genomeFilename = getGenomeId() + Globals.GENOME_FILE_EXTENSION;
        }
        return genomeFilename;
    }


    protected File showGenomeArchiveDirectoryChooser() {

        File directory = PreferenceManager.getInstance().getLastGenomeImportDirectory();
        File archiveName = new File(getGenomeId() + Globals.GENOME_FILE_EXTENSION);
        File file = FileDialogUtils.chooseFile("Save Genome File", directory, archiveName, FileDialogUtils.SAVE);

        if (file != null) {

            genomeFilename = file.getName();
            if (genomeFilename != null) {
                if (!genomeFilename.endsWith(Globals.GENOME_FILE_EXTENSION)) {
                    genomeFilename += Globals.GENOME_FILE_EXTENSION;
                    file = new File(file.getParentFile(), genomeFilename);
                }
                genomeArchiveLocation = file.getParentFile().getAbsolutePath();
            }
        }
        return file;
    }

    public boolean validateSelection() {

        try {

            if (!isIdValid()) {
                return false;
            }

            if (!isGenomeDisplayNameValid()) {
                genomeDisplayNameTextField.setText(null);
                return false;
            }

            if (!isFASTAFileValid()) {
                JOptionPane.showMessageDialog(this, "A fasta file is required!");
                return false;
            }

        } catch (Exception e) {
            logger.error("Error during Genome Builder validation!", e);
            return false;
        }
        return true;
    }

    private boolean isFASTAFileValid() {

        String file = fastaFileTextField.getText();
        return file != null && file.trim().length() > 0;
    }

    private Boolean isIdValid() {

        String id = getGenomeId();
        if (id == null || id.trim().length() == 0) {
            JOptionPane.showMessageDialog(this, "Genome ID is required");
            return false;
        }

        Collection<String> inUseIds = igv.getGenomeIds();
        if (inUseIds.contains(id)) {
            JOptionPane.showMessageDialog(this,
                    "The genome ID '" + id + "' is already in use - please select another!");
            return false;
        }
        return true;
    }

    private boolean isGenomeDisplayNameValid() {

        String displayName = getGenomeDisplayName();
        if (displayName == null || displayName.trim().length() == 0) {
            JOptionPane.showMessageDialog(this, "Genome name is required");
            return false;
        }

        Collection<String> inUseDisplayNames = igv.getGenomeDisplayNames();

        if (inUseDisplayNames.contains(displayName)) {
            JOptionPane.showMessageDialog(this,
                    "The genome name '" + displayName + "' is already in use - please select another!");
            return false;
        }
        return true;
    }

    /**
     * This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    // Generated using JFormDesigner non-commercial license
    private void initComponents() {
        jPanel1 = new JPanel();
        vSpacer3 = new JPanel(null);
        panel2 = new JPanel();
        genomeDisplayNameLabel2 = new JLabel();
        idField = new JTextField();
        genomeDisplayNameLabel = new JLabel();
        genomeDisplayNameTextField = new JTextField();
        fastaFileLabel = new JLabel();
        fastaFileTextField = new JTextField();
        fastaFileButton = new JButton();
        vSpacer1 = new JPanel(null);
        panel1 = new JPanel();
        cytobandFileLabel = new JLabel();
        cytobandFileTextField = new JTextField();
        cytobandFileButton = new JButton();
        refFlatFileLabel = new JLabel();
        refFlatFileTextField = new JTextField();
        refFlatFileButton = new JButton();
        refFlatFileLabel2 = new JLabel();
        chrAliasField = new JTextField();
        chrAliasButton = new JButton();
        vSpacer2 = new JPanel(null);
        label2 = new JLabel();
        label3 = new JLabel();
        sequenceURLField = new JTextField();

        //======== this ========
        setFont(new Font("Tahoma", Font.ITALIC, 12));
        setMaximumSize(new Dimension(900, 500));
        setMinimumSize(new Dimension(400, 300));
        setPreferredSize(new Dimension(700, 400));
        setLayout(new BorderLayout());

        //======== jPanel1 ========
        {
            jPanel1.setBorder(null);
            jPanel1.setLayout(new GridBagLayout());
            ((GridBagLayout)jPanel1.getLayout()).columnWidths = new int[] {15, 0, 0, 578, 83, 0};
            ((GridBagLayout)jPanel1.getLayout()).rowHeights = new int[] {15, 0, 25, 0, 0, 25, 0};
            ((GridBagLayout)jPanel1.getLayout()).columnWeights = new double[] {0.0, 0.0, 0.0, 0.0, 0.0, 1.0E-4};
            jPanel1.add(vSpacer3, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
                GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                new Insets(0, 0, 5, 5), 0, 0));

            //======== panel2 ========
            {
                panel2.setBorder(new TitledBorder(null, "Required", TitledBorder.LEADING, TitledBorder.DEFAULT_POSITION,
                    new Font("Lucida Grande", Font.BOLD, 13)));
                panel2.setLayout(new GridBagLayout());
                ((GridBagLayout)panel2.getLayout()).columnWidths = new int[] {0, 0, 578, 83, 0};
                ((GridBagLayout)panel2.getLayout()).columnWeights = new double[] {0.0, 0.0, 0.0, 0.0, 1.0E-4};

                //---- genomeDisplayNameLabel2 ----
                genomeDisplayNameLabel2.setText("ID");
                genomeDisplayNameLabel2.setToolTipText("Unique identifier (e.g. hg18)");
                genomeDisplayNameLabel2.setMaximumSize(new Dimension(84, 16));
                genomeDisplayNameLabel2.setMinimumSize(new Dimension(84, 16));
                genomeDisplayNameLabel2.setPreferredSize(new Dimension(84, 16));
                panel2.add(genomeDisplayNameLabel2, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 5, 5), 0, 0));

                //---- idField ----
                idField.setToolTipText("A uniqe identifier for the genome");
                panel2.add(idField, new GridBagConstraints(1, 0, 2, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 5, 5), 0, 0));

                //---- genomeDisplayNameLabel ----
                genomeDisplayNameLabel.setText("Name");
                panel2.add(genomeDisplayNameLabel, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 5, 5), 0, 0));

                //---- genomeDisplayNameTextField ----
                genomeDisplayNameTextField.setToolTipText("The user-readable name of the genome");
                genomeDisplayNameTextField.setPreferredSize(new Dimension(400, 28));
                genomeDisplayNameTextField.setMinimumSize(new Dimension(25, 28));
                panel2.add(genomeDisplayNameTextField, new GridBagConstraints(1, 1, 2, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 5, 5), 0, 0));

                //---- fastaFileLabel ----
                fastaFileLabel.setText("Fasta file");
                panel2.add(fastaFileLabel, new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 0, 5), 0, 0));

                //---- fastaFileTextField ----
                fastaFileTextField.setToolTipText("A FASTA data file");
                fastaFileTextField.setPreferredSize(new Dimension(400, 28));
                fastaFileTextField.setMinimumSize(new Dimension(25, 28));
                panel2.add(fastaFileTextField, new GridBagConstraints(1, 2, 2, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 0, 5), 0, 0));

                //---- fastaFileButton ----
                fastaFileButton.setLabel("...");
                fastaFileButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        fastaFileButtonActionPerformed(e);
                    }
                });
                panel2.add(fastaFileButton, new GridBagConstraints(3, 2, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 0, 0), 0, 0));
            }
            jPanel1.add(panel2, new GridBagConstraints(1, 1, 4, 1, 0.0, 0.0,
                GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                new Insets(0, 0, 5, 0), 0, 0));
            jPanel1.add(vSpacer1, new GridBagConstraints(1, 2, 1, 1, 0.0, 0.0,
                GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                new Insets(0, 0, 5, 5), 0, 0));

            //======== panel1 ========
            {
                panel1.setBorder(new TitledBorder(null, "Optional", TitledBorder.LEADING, TitledBorder.DEFAULT_POSITION,
                    new Font("Lucida Grande", Font.BOLD, 13)));
                panel1.setLayout(new GridBagLayout());
                ((GridBagLayout)panel1.getLayout()).columnWidths = new int[] {0, 0, 578, 83, 0};
                ((GridBagLayout)panel1.getLayout()).rowHeights = new int[] {0, 0, 0, 25, 0, 0};
                ((GridBagLayout)panel1.getLayout()).columnWeights = new double[] {0.0, 0.0, 0.0, 0.0, 1.0E-4};

                //---- cytobandFileLabel ----
                cytobandFileLabel.setText("Cytoband file");
                panel1.add(cytobandFileLabel, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 5, 5), 0, 0));

                //---- cytobandFileTextField ----
                cytobandFileTextField.setToolTipText("A cytoband data file");
                cytobandFileTextField.setPreferredSize(new Dimension(400, 28));
                cytobandFileTextField.setMinimumSize(new Dimension(25, 28));
                panel1.add(cytobandFileTextField, new GridBagConstraints(1, 0, 2, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 5, 5), 0, 0));

                //---- cytobandFileButton ----
                cytobandFileButton.setLabel("...");
                cytobandFileButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        cytobandFileButtonActionPerformed(e);
                    }
                });
                panel1.add(cytobandFileButton, new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 5, 0), 0, 0));

                //---- refFlatFileLabel ----
                refFlatFileLabel.setText("Gene file");
                panel1.add(refFlatFileLabel, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 5, 5), 0, 0));

                //---- refFlatFileTextField ----
                refFlatFileTextField.setToolTipText("An annotation file");
                refFlatFileTextField.setPreferredSize(new Dimension(400, 28));
                refFlatFileTextField.setMinimumSize(new Dimension(25, 28));
                panel1.add(refFlatFileTextField, new GridBagConstraints(1, 1, 2, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 5, 5), 0, 0));

                //---- refFlatFileButton ----
                refFlatFileButton.setLabel("...");
                refFlatFileButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        refFlatFileButtonActionPerformed(e);
                    }
                });
                panel1.add(refFlatFileButton, new GridBagConstraints(3, 1, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 5, 0), 0, 0));

                //---- refFlatFileLabel2 ----
                refFlatFileLabel2.setText("Alias file");
                panel1.add(refFlatFileLabel2, new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 5, 5), 0, 0));

                //---- chrAliasField ----
                chrAliasField.setPreferredSize(new Dimension(400, 28));
                chrAliasField.setMinimumSize(new Dimension(25, 28));
                panel1.add(chrAliasField, new GridBagConstraints(1, 2, 2, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 5, 5), 0, 0));

                //---- chrAliasButton ----
                chrAliasButton.setLabel("...");
                chrAliasButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        chrAliasButtonActionPerformed(e);
                    }
                });
                panel1.add(chrAliasButton, new GridBagConstraints(3, 2, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 5, 0), 0, 0));
                panel1.add(vSpacer2, new GridBagConstraints(0, 3, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 5, 5), 0, 0));

                //---- label2 ----
                label2.setText("Supply a fasta URL if defining a web-hosted genome (optional, not common).");
                label2.setFont(new Font("Lucida Grande", Font.ITALIC, 13));
                panel1.add(label2, new GridBagConstraints(0, 4, 4, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 5, 0), 0, 0));

                //---- label3 ----
                label3.setText("Fasta URL");
                panel1.add(label3, new GridBagConstraints(0, 5, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 0, 5), 0, 0));

                //---- sequenceURLField ----
                sequenceURLField.setToolTipText("A refFlat gene file");
                sequenceURLField.setPreferredSize(new Dimension(400, 28));
                sequenceURLField.setMinimumSize(new Dimension(25, 28));
                panel1.add(sequenceURLField, new GridBagConstraints(1, 5, 2, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 0, 5), 0, 0));
            }
            jPanel1.add(panel1, new GridBagConstraints(1, 3, 4, 1, 0.0, 0.0,
                GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                new Insets(0, 0, 5, 0), 0, 0));
        }
        add(jPanel1, BorderLayout.CENTER);
    }// </editor-fold>//GEN-END:initComponents

    private void chrAliasButtonActionPerformed(ActionEvent e) {
        File directory = PreferenceManager.getInstance().getDefineGenomeInputDirectory();
        File file = FileDialogUtils.chooseFile("Select Chromosome Alias File", directory, FileDialogUtils.LOAD);
        if (file != null) {
            chrAliasField.setText(file.getAbsolutePath());
            PreferenceManager.getInstance().setDefineGenomeInputDirectory(file.getParentFile());
        }
    }


    private void cytobandFileButtonActionPerformed(java.awt.event.ActionEvent evt) {
        File directory = PreferenceManager.getInstance().getDefineGenomeInputDirectory();
        File file = FileDialogUtils.chooseFile("Select Cytoband File", directory, FileDialogUtils.LOAD);
        if (file != null) {
            cytobandFileTextField.setText(file.getAbsolutePath());
            PreferenceManager.getInstance().setDefineGenomeInputDirectory(file.getParentFile());
        }
    }

    private void refFlatFileButtonActionPerformed(java.awt.event.ActionEvent evt) {
        File directory = PreferenceManager.getInstance().getDefineGenomeInputDirectory();
        File file = FileDialogUtils.chooseFile("Select Annotation File", directory, FileDialogUtils.LOAD);
        if (file != null) {
            refFlatFileTextField.setText(file.getAbsolutePath());
            PreferenceManager.getInstance().setDefineGenomeInputDirectory(file.getParentFile());
        }
    }

    private void fastaFileButtonActionPerformed(java.awt.event.ActionEvent evt) {
        File directory = PreferenceManager.getInstance().getDefineGenomeInputDirectory();
        File file = FileDialogUtils.chooseFile("Select Fasta File", directory, FileDialogUtils.LOAD);
        if (file != null) {
            fastaFileTextField.setText(file.getAbsolutePath());
            PreferenceManager.getInstance().setDefineGenomeInputDirectory(file.getParentFile());
        }
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel jPanel1;
    private JPanel vSpacer3;
    private JPanel panel2;
    private JLabel genomeDisplayNameLabel2;
    private JTextField idField;
    private JLabel genomeDisplayNameLabel;
    private JTextField genomeDisplayNameTextField;
    private JLabel fastaFileLabel;
    private JTextField fastaFileTextField;
    private JButton fastaFileButton;
    private JPanel vSpacer1;
    private JPanel panel1;
    private JLabel cytobandFileLabel;
    private JTextField cytobandFileTextField;
    private JButton cytobandFileButton;
    private JLabel refFlatFileLabel;
    private JTextField refFlatFileTextField;
    private JButton refFlatFileButton;
    private JLabel refFlatFileLabel2;
    private JTextField chrAliasField;
    private JButton chrAliasButton;
    private JPanel vSpacer2;
    private JLabel label2;
    private JLabel label3;
    private JTextField sequenceURLField;
    // End of variables declaration//GEN-END:variables

}
