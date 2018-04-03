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


package org.broad.igv.ui.util;

import org.apache.log4j.Logger;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.sam.reader.AlignmentIndexer;
import org.broad.igv.sam.reader.FeatureIndex;
import org.broad.igv.tools.IgvTools;
import org.broad.igv.util.ResourceLocator;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;

import java.awt.*;
import java.awt.event.*;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import javax.swing.*;
import javax.swing.border.*;


public class IndexCreatorDialog extends JDialog {

    private static Logger log = Logger.getLogger(IndexCreatorDialog.class);

    File file;
    File idxFile;
    IndexWorker worker;
    FileType fileType;
    String introText;

    private enum FileType {
        SAM,
        TRIBBLE
    }


    public static IndexCreatorDialog createShowDialog(Frame parent, File baseFile, File newIdxFile, String introText) {
        final IndexCreatorDialog dialog = new IndexCreatorDialog(parent, true, baseFile, newIdxFile, introText);
        dialog.setLocationRelativeTo(parent);

        UIUtilities.invokeAndWaitOnEventThread(new Runnable() {
            @Override
            public void run() {
                dialog.setVisible(true);
            }
        });
        return dialog;
    }

    /**
     * Creates new form IndexCreatorDialog
     */
    public IndexCreatorDialog(java.awt.Frame parent, boolean modal, File file, File idxFile, String introText) {
        super(parent, modal);
        initComponents();
        this.introText = introText;

        jLabel1.setVisible(false);

        this.file = file;
        this.idxFile = idxFile;
        this.determineFileType(file);
        if (this.fileType == null) {
            log.error("Cannot determine file type for " + file.getAbsolutePath());
        }

        String txt = introText;

        this.introTextArea.setText(txt);

        this.introTextArea.setBorder(BorderFactory.createEmptyBorder());

        switch (this.fileType) {
            case SAM:
                worker = new SamIndexWorker();
                break;
            case TRIBBLE:
                worker = new TribbleIndexWorker();
                break;
        }
    }

    private void determineFileType(File file) {

        String filename = file.getName();
        if (filename.toLowerCase().endsWith(".sam")) {
            fileType = FileType.SAM;
        } else if (CodecFactory.hasCodec(new ResourceLocator(file.getPath()), null)) {
            fileType = FileType.TRIBBLE;
        } else {
            fileType = null;
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
                if (fileType == FileType.TRIBBLE) {
                    IndexCreatorDialog.this.progressBar.setIndeterminate(true);
                    jLabel1.setText("Creating index...");
                }
            }
        }
    }

    private void cancelButtonActionPerformed(java.awt.event.ActionEvent evt) {

        if (worker != null && worker.isStarted && !(worker.isDone() || worker.isCancelled())) {
            worker.isStarted = false;
            worker.cancel(true);  // <+ doing this before execution starts will raise an error
        }
        setVisible(false);
    }

    public class SamIndexWorker extends IndexWorker<FeatureIndex> {
        @Override
        protected FeatureIndex doInBackground() throws Exception {
            AlignmentIndexer indexer = AlignmentIndexer.getInstance(file, progressBar, this);
            return indexer.createSamIndex(idxFile, 16000);
        }
    }

    private class TribbleIndexWorker extends IndexWorker<Index> {
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
    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        scrollPane1 = new JScrollPane();
        introTextArea = new JTextPane();
        timeRemainingLabel = new JLabel();
        jLabel1 = new JLabel();
        progressBar = new JProgressBar();
        buttonBar = new JPanel();
        goButton = new JButton();
        cancelButton = new JButton();

        //======== this ========
        setResizable(false);
        setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(null);

                //======== scrollPane1 ========
                {
                    scrollPane1.setViewportView(introTextArea);
                }
                contentPanel.add(scrollPane1);
                scrollPane1.setBounds(0, 0, 465, 260);
                contentPanel.add(timeRemainingLabel);
                timeRemainingLabel.setBounds(255, 265, 210, 31);

                //---- jLabel1 ----
                jLabel1.setHorizontalTextPosition(SwingConstants.RIGHT);
                contentPanel.add(jLabel1);
                jLabel1.setBounds(0, 265, 210, 31);
                contentPanel.add(progressBar);
                progressBar.setBounds(0, 305, 465, progressBar.getPreferredSize().height);

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
                buttonBar.setLayout(new GridBagLayout());
                ((GridBagLayout)buttonBar.getLayout()).columnWidths = new int[] {0, 85, 80};
                ((GridBagLayout)buttonBar.getLayout()).columnWeights = new double[] {1.0, 0.0, 0.0};

                //---- goButton ----
                goButton.setText("Go");
                goButton.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        goButtonActionPerformed(e);
                    }
                });
                buttonBar.add(goButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
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
        setSize(500, 415);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JScrollPane scrollPane1;
    private JTextPane introTextArea;
    private JLabel timeRemainingLabel;
    private JLabel jLabel1;
    private JProgressBar progressBar;
    private JPanel buttonBar;
    private JButton goButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
