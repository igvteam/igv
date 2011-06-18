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
        cytobandFileTextField = new JTextField();
        cytobandFileButton = new JButton();
        refFlatFileLabel = new JLabel();
        refFlatFileButton = new JButton();
        genomeDisplayNameTextField = new JTextField();
        fastaFileLabel = new JLabel();
        fastaFileButton = new JButton();
        cytobandFileLabel = new JLabel();
        fastaFileTextField = new JTextField();
        refFlatFileTextField = new JTextField();
        genomeDisplayNameLabel = new JLabel();
        refFlatFileLabel2 = new JLabel();
        chrAliasField = new JTextField();
        chrAliasButton = new JButton();
        idField = new JTextField();
        genomeDisplayNameLabel2 = new JLabel();
        label1 = new JLabel();
        sequenceDirectoryCB = new JCheckBox();
        label2 = new JLabel();
        label3 = new JLabel();
        sequenceURLField = new JTextField();
        jLabel1 = new JLabel();

        //======== this ========
        setFont(new Font("Tahoma", Font.ITALIC, 12));
        setMaximumSize(new Dimension(900, 500));
        setMinimumSize(new Dimension(400, 300));
        setPreferredSize(new Dimension(700, 400));
        setLayout(null);

        //======== jPanel1 ========
        {
            jPanel1.setBorder(null);
            jPanel1.setLayout(null);

            //---- cytobandFileTextField ----
            cytobandFileTextField.setToolTipText("A cytoband data file");
            cytobandFileTextField.setPreferredSize(new Dimension(400, 28));
            cytobandFileTextField.setMinimumSize(new Dimension(25, 28));
            jPanel1.add(cytobandFileTextField);
            cytobandFileTextField.setBounds(105, 135, 608, 29);

            //---- cytobandFileButton ----
            cytobandFileButton.setLabel("...");
            cytobandFileButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    cytobandFileButtonActionPerformed(e);
                }
            });
            jPanel1.add(cytobandFileButton);
            cytobandFileButton.setBounds(725, 130, 50, cytobandFileButton.getPreferredSize().height);

            //---- refFlatFileLabel ----
            refFlatFileLabel.setText("Gene file");
            jPanel1.add(refFlatFileLabel);
            refFlatFileLabel.setBounds(15, 170, 87, 29);

            //---- refFlatFileButton ----
            refFlatFileButton.setLabel("...");
            refFlatFileButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    refFlatFileButtonActionPerformed(e);
                }
            });
            jPanel1.add(refFlatFileButton);
            refFlatFileButton.setBounds(725, 165, 50, refFlatFileButton.getPreferredSize().height);

            //---- genomeDisplayNameTextField ----
            genomeDisplayNameTextField.setToolTipText("The user-readable name of the genome");
            genomeDisplayNameTextField.setPreferredSize(new Dimension(400, 28));
            genomeDisplayNameTextField.setMinimumSize(new Dimension(25, 28));
            jPanel1.add(genomeDisplayNameTextField);
            genomeDisplayNameTextField.setBounds(105, 40, 288, genomeDisplayNameTextField.getPreferredSize().height);

            //---- fastaFileLabel ----
            fastaFileLabel.setText("Fasta file *");
            jPanel1.add(fastaFileLabel);
            fastaFileLabel.setBounds(15, 100, 87, 29);

            //---- fastaFileButton ----
            fastaFileButton.setLabel("...");
            fastaFileButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    fastaFileButtonActionPerformed(e);
                }
            });
            jPanel1.add(fastaFileButton);
            fastaFileButton.setBounds(725, 95, 50, fastaFileButton.getPreferredSize().height);

            //---- cytobandFileLabel ----
            cytobandFileLabel.setText("Cytoband file");
            jPanel1.add(cytobandFileLabel);
            cytobandFileLabel.setBounds(15, 135, 87, 29);

            //---- fastaFileTextField ----
            fastaFileTextField.setToolTipText("A FASTA data file");
            fastaFileTextField.setPreferredSize(new Dimension(400, 28));
            fastaFileTextField.setMinimumSize(new Dimension(25, 28));
            jPanel1.add(fastaFileTextField);
            fastaFileTextField.setBounds(105, 100, 608, 29);

            //---- refFlatFileTextField ----
            refFlatFileTextField.setToolTipText("An annotation file");
            refFlatFileTextField.setPreferredSize(new Dimension(400, 28));
            refFlatFileTextField.setMinimumSize(new Dimension(25, 28));
            jPanel1.add(refFlatFileTextField);
            refFlatFileTextField.setBounds(105, 170, 608, 29);

            //---- genomeDisplayNameLabel ----
            genomeDisplayNameLabel.setText("Name *");
            jPanel1.add(genomeDisplayNameLabel);
            genomeDisplayNameLabel.setBounds(15, 40, 87, 28);

            //---- refFlatFileLabel2 ----
            refFlatFileLabel2.setText("Alias file");
            jPanel1.add(refFlatFileLabel2);
            refFlatFileLabel2.setBounds(15, 220, 87, 29);

            //---- chrAliasField ----
            chrAliasField.setPreferredSize(new Dimension(400, 28));
            chrAliasField.setMinimumSize(new Dimension(25, 28));
            jPanel1.add(chrAliasField);
            chrAliasField.setBounds(105, 220, 608, 29);

            //---- chrAliasButton ----
            chrAliasButton.setLabel("...");
            chrAliasButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    chrAliasButtonActionPerformed(e);
                }
            });
            jPanel1.add(chrAliasButton);
            chrAliasButton.setBounds(725, 220, 50, chrAliasButton.getPreferredSize().height);

            //---- idField ----
            idField.setToolTipText("A uniqe identifier for the genome");
            jPanel1.add(idField);
            idField.setBounds(105, 5, 148, idField.getPreferredSize().height);

            //---- genomeDisplayNameLabel2 ----
            genomeDisplayNameLabel2.setText("ID *");
            jPanel1.add(genomeDisplayNameLabel2);
            genomeDisplayNameLabel2.setBounds(15, 5, 87, 28);

            //---- label1 ----
            label1.setText("(unique id, e.g. hg18)");
            jPanel1.add(label1);
            label1.setBounds(270, 5, 166, 28);

            //---- sequenceDirectoryCB ----
            sequenceDirectoryCB.setText("Fasta file is a directory");
            sequenceDirectoryCB.setFont(new Font("Lucida Grande", Font.ITALIC, 13));
            jPanel1.add(sequenceDirectoryCB);
            sequenceDirectoryCB.setBounds(105, 75, 546, sequenceDirectoryCB.getPreferredSize().height);

            //---- label2 ----
            label2.setText("Supply a sequence URL if defining a web-hosted genome (optional, not common).  See user guide for more details. ");
            label2.setFont(new Font("Lucida Grande", Font.ITALIC, 13));
            jPanel1.add(label2);
            label2.setBounds(new Rectangle(new Point(15, 275), label2.getPreferredSize()));

            //---- label3 ----
            label3.setText("Sequence URL");
            jPanel1.add(label3);
            label3.setBounds(15, 295, label3.getPreferredSize().width, 28);

            //---- sequenceURLField ----
            sequenceURLField.setToolTipText("A refFlat gene file");
            sequenceURLField.setPreferredSize(new Dimension(400, 28));
            sequenceURLField.setMinimumSize(new Dimension(25, 28));
            jPanel1.add(sequenceURLField);
            sequenceURLField.setBounds(110, 295, 608, sequenceURLField.getPreferredSize().height);

            { // compute preferred size
                Dimension preferredSize = new Dimension();
                for(int i = 0; i < jPanel1.getComponentCount(); i++) {
                    Rectangle bounds = jPanel1.getComponent(i).getBounds();
                    preferredSize.width = Math.max(bounds.x + bounds.width, preferredSize.width);
                    preferredSize.height = Math.max(bounds.y + bounds.height, preferredSize.height);
                }
                Insets insets = jPanel1.getInsets();
                preferredSize.width += insets.right;
                preferredSize.height += insets.bottom;
                jPanel1.setMinimumSize(preferredSize);
                jPanel1.setPreferredSize(preferredSize);
            }
        }
        add(jPanel1);
        jPanel1.setBounds(10, 0, 810, jPanel1.getPreferredSize().height);

        //---- jLabel1 ----
        jLabel1.setFont(new Font("Lucida Sans Unicode", Font.ITALIC, 12));
        jLabel1.setText("<html>* required.  <br>The sequence file (required) can be a FASTA file, a directory of FASTA files, or a zip of FASTA files. Optionally, specify a cytoband file to display the chromosome ideogram and an annotation file to display the gene track. See the documentation for descriptions of supported annotation formats.");
        jLabel1.setVerticalAlignment(SwingConstants.TOP);
        jLabel1.setBorder(new EmptyBorder(5, 25, 25, 5));
        add(jLabel1);
        jLabel1.setBounds(5, 340, 785, 150);

        { // compute preferred size
            Dimension preferredSize = new Dimension();
            for(int i = 0; i < getComponentCount(); i++) {
                Rectangle bounds = getComponent(i).getBounds();
                preferredSize.width = Math.max(bounds.x + bounds.width, preferredSize.width);
                preferredSize.height = Math.max(bounds.y + bounds.height, preferredSize.height);
            }
            Insets insets = getInsets();
            preferredSize.width += insets.right;
            preferredSize.height += insets.bottom;
            setMinimumSize(preferredSize);
            setPreferredSize(preferredSize);
        }
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
        boolean chooseDir = sequenceDirectoryCB.isSelected();
        File directory = PreferenceManager.getInstance().getDefineGenomeInputDirectory();
        File file = chooseDir ? FileDialogUtils.chooseDirectory("Select Fasta Directory", directory) :
                FileDialogUtils.chooseFile("Select Fasta File", directory, FileDialogUtils.LOAD);
        if (file != null) {
            fastaFileTextField.setText(file.getAbsolutePath());
            PreferenceManager.getInstance().setDefineGenomeInputDirectory(file.getParentFile());
        }
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel jPanel1;
    private JTextField cytobandFileTextField;
    private JButton cytobandFileButton;
    private JLabel refFlatFileLabel;
    private JButton refFlatFileButton;
    private JTextField genomeDisplayNameTextField;
    private JLabel fastaFileLabel;
    private JButton fastaFileButton;
    private JLabel cytobandFileLabel;
    private JTextField fastaFileTextField;
    private JTextField refFlatFileTextField;
    private JLabel genomeDisplayNameLabel;
    private JLabel refFlatFileLabel2;
    private JTextField chrAliasField;
    private JButton chrAliasButton;
    private JTextField idField;
    private JLabel genomeDisplayNameLabel2;
    private JLabel label1;
    private JCheckBox sequenceDirectoryCB;
    private JLabel label2;
    private JLabel label3;
    private JTextField sequenceURLField;
    private JLabel jLabel1;
    // End of variables declaration//GEN-END:variables

}
