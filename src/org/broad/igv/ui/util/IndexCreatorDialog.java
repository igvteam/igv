/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * SamIndexCreatorDialog.java
 *
 * Created on Feb 11, 2009, 1:45:19 PM
 */
package org.broad.igv.ui.util;

import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.sam.reader.AlignmentIndexer;
import org.broad.igv.sam.reader.FeatureIndex;
import org.broad.igv.tools.IgvTools;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.TribbleException;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;

import javax.swing.*;
import java.awt.*;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;

/**
 * Dialog for asking the user if they want to create an index, and
 * displaying progress if they do.
 *
 * @author jacob
 */

public class IndexCreatorDialog extends JDialog {

    File file;
    File idxFile;
    IndexWorker worker;

    FileType fileType;

    private enum FileType {
        SAM,
        VCF
    }

    static String introText = "An index file for @filename could not " +
            "be located. An index is recommended to view @filetype files in IGV.  " +
            "Click \"Go\" to create one now.";

    public static IndexCreatorDialog createShowDialog(Frame parent, File baseFile, File newIdxFile) {
        IndexCreatorDialog dialog = new IndexCreatorDialog(parent, true, baseFile, newIdxFile);
        dialog.setLocationRelativeTo(parent);
        dialog.setVisible(true);
        return dialog;
    }

    /**
     * Creates new form IndexCreatorDialog
     */
    public IndexCreatorDialog(java.awt.Frame parent, boolean modal,
                              File file,
                              File idxFile) {
        super(parent, modal);
        initComponents();
        jLabel1.setVisible(false);

        this.file = file;
        this.idxFile = idxFile;
        this.determineFileType(file);
        if (this.fileType == null)
            throw new IllegalArgumentException("Cannot determine file type for " + file.getAbsolutePath());

        int timeEst = 1 + (int) Math.ceil(file.length() / 1000000000.0);

        String txt = introText.replace("@filename", file.getName()).replace(
                "@time", String.valueOf(timeEst)).replace("@filetype", this.fileType.name());

        this.introTextPane.setText(txt);

        this.introTextPane.setBorder(BorderFactory.createEmptyBorder());

        switch (this.fileType) {
            case SAM:
                worker = new SamIndexWorker();
                break;
            case VCF:
                worker = new VCFIndexWorker();
                break;
        }
    }

    private void determineFileType(File file) {
        for (FileType ft : FileType.values()) {
            if (file.getName().toLowerCase().endsWith(ft.name().toLowerCase())) {
                this.fileType = ft;
            }
        }
    }

    public Object getIndex() {
        if (worker == null || !worker.isDone()) {
            return null;
        } else {
            try {
                return worker.get();
            } catch (Exception ex) {
                MessageUtils.showMessage(ex.getMessage());
            }
            return null;
        }
    }

    public class SamIndexWorker extends IndexWorker<FeatureIndex> {
        @Override
        protected FeatureIndex doInBackground() throws Exception {
            AlignmentIndexer indexer = AlignmentIndexer.getInstance(file, progressBar, this);
            return indexer.createSamIndex(idxFile, 16000);
        }
    }

    private class VCFIndexWorker extends IndexWorker<Index> {
        @Override
        protected Index doInBackground() throws Exception {
            int binSize = IgvTools.LINEAR_BIN_SIZE;
            FeatureCodec codec = CodecFactory.getCodec(file.getAbsolutePath(), GenomeManager.getInstance().getCurrentGenome());
            if (codec != null) {
                try {
                    Index index = IndexFactory.createLinearIndex(file, codec, binSize);
                    if (index != null) {
                        IgvTools.writeTribbleIndex(index, idxFile.getAbsolutePath());
                    }
                    return index;
                } catch (TribbleException.MalformedFeatureFile e) {
                    StringBuffer buf = new StringBuffer();
                    buf.append("<html>Files must be sorted by start position prior to indexing.<br>");
                    buf.append(e.getMessage());
                    buf.append("<br><br>Note: igvtools can be used to sort the file, select \"File > Run igvtools...\".");
                    MessageUtils.showMessage(buf.toString());
                }
            } else {
                throw new DataLoadException("Unknown File Type", file.getAbsolutePath());
            }
            return null;
        }
    }

    public abstract class IndexWorker<I> extends SwingWorker<I, Void> {

        private boolean isStarted = false;

        @Override
        protected void done() {
            setVisible(false);
        }

        public void setTimeRemaining(long timeInMillis) {
            final int timeRemaining = (int) (timeInMillis / (60 * 1000));
            UIUtilities.invokeOnEventThread(new Runnable() {

                public void run() {
                    String txt = String.valueOf(timeRemaining) + " minutes";
                    if (timeRemaining == 1) {
                        txt = "1 minute";
                    } else if (timeRemaining < 1) {
                        txt = " < 1 minute";
                    }
                    timeRemainingLabel.setText(txt);
                }
            });

        }


    }

    /**
     * ProgressListener listens to "progress" property
     * changes in the SwingWorkers that search and load
     * images.
     */
    class ProgressListener implements PropertyChangeListener {
        // prevent creation without providing a progress bar
        private ProgressListener() {
        }

        ProgressListener(JProgressBar progressBar) {
            this.progressBar = progressBar;
            this.progressBar.setValue(0);
        }

        public void propertyChange(PropertyChangeEvent evt) {
            String strPropertyName = evt.getPropertyName();
            if ("progress".equals(strPropertyName)) {
                progressBar.setIndeterminate(false);
                int progress = (Integer) evt.getNewValue();
                progressBar.setValue(progress);
            }
        }

        private JProgressBar progressBar;
    }

    /**
     * This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jLabel1 = new JLabel();
        progressBar = new JProgressBar();
        timeRemainingLabel = new JLabel();
        goButton = new JButton();
        cancelButton = new JButton();
        jScrollPane2 = new JScrollPane();
        introTextPane = new JTextPane();

        setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);

        jLabel1.setText("Estimated time remaining: ");

        goButton.setText("Go");
        goButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                goButtonActionPerformed(evt);
            }
        });

        cancelButton.setText("Cancel");
        cancelButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                cancelButtonActionPerformed(evt);
            }
        });

        introTextPane.setBackground(getParent().getBackground());
        introTextPane.setBorder(BorderFactory.createEmptyBorder(1, 1, 1, 1));
        introTextPane.setEditable(false);
        introTextPane.setText("An index file for [filename goes here] could not be located.  An index is required to view alignments in IGV.  Click \"Go\" to create one now.  This will take approximately [time goes here] to complete.");
        introTextPane.setFocusable(false);
        jScrollPane2.setViewportView(introTextPane);

        org.jdesktop.layout.GroupLayout layout = new org.jdesktop.layout.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
                layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                        .add(layout.createSequentialGroup()
                                .add(30, 30, 30)
                                .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                                        .add(jScrollPane2, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 343, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                                        .add(layout.createSequentialGroup()
                                                .add(jLabel1)
                                                .add(18, 18, 18)
                                                .add(timeRemainingLabel, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 103, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                                        .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.TRAILING)
                                                .add(layout.createSequentialGroup()
                                                        .add(cancelButton)
                                                        .add(18, 18, 18)
                                                        .add(goButton))
                                                .add(progressBar, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 351, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)))
                                .addContainerGap(26, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
                layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                        .add(org.jdesktop.layout.GroupLayout.TRAILING, layout.createSequentialGroup()
                                .addContainerGap(23, Short.MAX_VALUE)
                                .add(jScrollPane2, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 132, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                                .add(35, 35, 35)
                                .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                                        .add(jLabel1, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 27, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                                        .add(timeRemainingLabel, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 16, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                                .add(progressBar, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                                .add(18, 18, 18)
                                .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                                        .add(goButton)
                                        .add(cancelButton))
                                .addContainerGap())
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void goButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_goButtonActionPerformed

        if (worker.isDone() || worker.isCancelled()) {
            setVisible(false);
        } else {
            if (!worker.isStarted) {
                goButton.setEnabled(false);
                worker.isStarted = true;
                worker.execute();
                jLabel1.setVisible(true);
                //Haven't worked out how to publish progress yet, just going to set it to indeterminate
                if (fileType == FileType.VCF) {
                    IndexCreatorDialog.this.progressBar.setIndeterminate(true);
                    jLabel1.setText("Creating index...");
                }
            }
        }
    }

    private void cancelButtonActionPerformed(java.awt.event.ActionEvent evt) {

        if(worker != null && worker.isStarted && !(worker.isDone() || worker.isCancelled())) {
            worker.isStarted = false;
            worker.cancel(true);  // <+ doing this before execution starts will raise an error
        }
        setVisible(false);
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        java.awt.EventQueue.invokeLater(new Runnable() {

            public void run() {
                File samFile = new File("/Users/jrobinso/IGV/Sam/30DWM.7.sam");
                File idxFile = new File("/Users/jrobinso/IGV/Sam/30DWM.7.sai");
                IndexCreatorDialog dialog = new IndexCreatorDialog(
                        new JFrame(), true,
                        samFile, idxFile);
                dialog.addWindowListener(new java.awt.event.WindowAdapter() {

                    @Override
                    public void windowClosing(java.awt.event.WindowEvent e) {
                        System.exit(0);
                    }
                });
                dialog.setVisible(true);
            }
        });
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private JButton cancelButton;
    private JButton goButton;
    private JTextPane introTextPane;
    private JLabel jLabel1;
    private JScrollPane jScrollPane2;
    private JProgressBar progressBar;
    private JLabel timeRemainingLabel;
    // End of variables declaration//GEN-END:variables


}
