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
 * GenomeSelectionDialog.java
 *
 * Created on November 8, 2007, 3:51 PM
 */

package org.broad.igv.ui.commandbar;

import org.broad.igv.DirectoryManager;
import org.broad.igv.event.GenomeResetEvent;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.feature.genome.GenomeDownloadUtils;
import org.broad.igv.feature.genome.GenomeListItem;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.LongRunningTask;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Dialog box for selecting genomes. User can choose from a list,
 * which is filtered according to text box input
 */
public class HostedGenomeSelectionDialog extends org.broad.igv.ui.IGVDialog {

    private static Logger log = LogManager.getLogger(HostedGenomeSelectionDialog.class);

    private JTextField genomeFilter;
    private JList<GenomeListItem> genomeList;
    private JCheckBox downloadSequenceCB;
    private JCheckBox downloadAnnotationsCB;
    private boolean isCanceled;
    private List<GenomeListItem> allListItems;
    private DefaultListModel genomeListModel;


    /**
     * Open a selection list to load a genome from the server.   This method is static because its used by multiple
     * UI elements  (menu bar and genome selection pulldown).
     */
    public static void selectHostedGenome() {

        Runnable showDialog = () -> {

            Collection<GenomeListItem> inputListItems = GenomeListManager.getInstance().getHostedGenomeList();
            if (inputListItems == null) {
                return;
            }
            HostedGenomeSelectionDialog dialog = new HostedGenomeSelectionDialog(IGV.getInstance().getMainFrame(), inputListItems);
            UIUtilities.invokeAndWaitOnEventThread(() -> dialog.setVisible(true));

            if (dialog.isCanceled()) {
                IGVEventBus.getInstance().post(new GenomeResetEvent());
            } else {

                GenomeListItem selectedItem = dialog.getSelectedValue();
                if (selectedItem == null) return;

                boolean downloadSequence = dialog.isDownloadSequence();
                boolean downloadAnnotations = dialog.isDownloadAnnotations();

                File downloadPath = null;
                if (downloadSequence || downloadAnnotations || selectedItem.getPath().endsWith(".genome")) {
                    downloadPath = GenomeManager.getInstance().downloadGenome(selectedItem, downloadSequence, downloadAnnotations);
                }


                try {
                    if (downloadPath != null) {
                        GenomeManager.getInstance().loadGenome(downloadPath.getAbsolutePath());
                    } else {
                        GenomeManager.getInstance().loadGenome(selectedItem.getPath());
                    }

                    // Legacy cleanup - json takes precedence over ".genome"
                    if (selectedItem.getPath().endsWith(".json")) {
                        removeDotGenomeFile(selectedItem.getId());
                    }
                } catch (IOException e) {
                    MessageUtils.showErrorMessage("Error loading genome " + selectedItem.getDisplayableName(), e);
                    log.error("Error loading genome " + selectedItem.getDisplayableName(), e);
                }
            }
        };

        if (SwingUtilities.isEventDispatchThread()) {
            LongRunningTask.submit(showDialog);
        } else {
            showDialog.run();
        }
    }

    private static void removeDotGenomeFile(String id) {
        try {
            File dotGenomeFile = new File(DirectoryManager.getGenomeCacheDirectory(), id + ".genome");
            if (dotGenomeFile.exists()) {
                dotGenomeFile.delete();
            }
        } catch (Exception e) {
            // If anything goes wrong, just log it, this cleanup is not essential
            log.error("Error deleting .genome file", e);
        }
    }

    /**
     * @param parent
     */
    private HostedGenomeSelectionDialog(java.awt.Frame parent, Collection<GenomeListItem> inputListItems) {
        super(parent);
        initComponents();
        initData(inputListItems);
    }

    private void initData(Collection<GenomeListItem> inputListItems) {
        this.allListItems = new ArrayList<>(inputListItems);
        String filterText = genomeFilter.getText().trim().toLowerCase();
        rebuildHostedGenomeList(filterText);
    }

    /**
     * Filter the list of displayed genomes with the text the user entered.
     */
    private void rebuildHostedGenomeList(String filterText) {
        if (genomeListModel == null) {
            genomeListModel = new DefaultListModel();
            UIUtilities.invokeOnEventThread(() -> genomeList.setModel(genomeListModel));
        }
        genomeListModel.clear();
        filterText = filterText.toLowerCase().trim();
        for (GenomeListItem listItem : allListItems) {
            if (listItem.getDisplayableName().toLowerCase().contains(filterText)) {
                genomeListModel.addElement(listItem);
            }
        }
    }


    private void genomeEntryKeyReleased(KeyEvent e) {
        rebuildHostedGenomeList(genomeFilter.getText());
    }

    public GenomeListItem getSelectedValue() {
        return isCanceled ? null : genomeList.getSelectedValue();
    }

    public boolean isDownloadSequence() {
        return downloadSequenceCB.isSelected();
    }

    public boolean isDownloadAnnotations() {
        return downloadAnnotationsCB.isSelected();
    }

    public boolean isCanceled() {
        return isCanceled;
    }

    private void cancelButtonActionPerformed(java.awt.event.ActionEvent evt) {
        isCanceled = true;
        setVisible(false);
        dispose();
    }

    private void okButtonActionPerformed(java.awt.event.ActionEvent evt) {
        isCanceled = false;
        setVisible(false);
        dispose();
    }


    private void configureDownloadButtons(GenomeListItem item) {
        if (item != null) {
            final boolean sequenceDownloadable = GenomeDownloadUtils.isSequenceDownloadable(item);
            final boolean annotationsDownloadable = GenomeDownloadUtils.isAnnotationsDownloadable(item);
            downloadSequenceCB.setEnabled(sequenceDownloadable);
            downloadAnnotationsCB.setEnabled(annotationsDownloadable);
        }
    }


    private void initComponents() {

        //======== this ========
        setModal(true);
        setTitle("Hosted Genomes");
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========

        JPanel dialogPane = new JPanel();
        dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
        dialogPane.setPreferredSize(new Dimension(350, 500));
        dialogPane.setLayout(new BorderLayout());

        //======== contentPanel ========
        JPanel contentPanel = new JPanel();
        contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.Y_AXIS));

        //---- textArea1 ----
        JTextArea textArea1 = new JTextArea();
        textArea1.setText("Selected genomes will be downloaded and added to the genome dropdown list.");
        textArea1.setLineWrap(true);
        textArea1.setWrapStyleWord(true);
        textArea1.setBackground(UIManager.getColor("Button.background"));
        textArea1.setRows(2);
        textArea1.setMaximumSize(new Dimension(2147483647, 60));
        textArea1.setRequestFocusEnabled(false);
        textArea1.setEditable(false);
        contentPanel.add(textArea1);

        //======== filterPanel ========
        JPanel filterPanel = new JPanel();
        filterPanel.setMaximumSize(new Dimension(2147483647, 28));
        filterPanel.setLayout(new GridBagLayout());
        ((GridBagLayout) filterPanel.getLayout()).columnWidths = new int[]{0, 0, 0};
        ((GridBagLayout) filterPanel.getLayout()).rowHeights = new int[]{0, 0};
        ((GridBagLayout) filterPanel.getLayout()).columnWeights = new double[]{1.0, 1.0, 1.0E-4};
        ((GridBagLayout) filterPanel.getLayout()).rowWeights = new double[]{1.0, 1.0E-4};

        //---- label1 ----
        JLabel filterLabel = new JLabel();
        filterLabel.setText("Filter:");
        filterLabel.setLabelFor(genomeFilter);
        filterLabel.setRequestFocusEnabled(false);
        filterPanel.add(filterLabel, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0,
                GridBagConstraints.WEST, GridBagConstraints.VERTICAL,
                new Insets(0, 0, 0, 0), 0, 0));

        //---- genomeFilter ----
        genomeFilter = new JTextField();
        genomeFilter.setToolTipText("Filter genome list");
        genomeFilter.setPreferredSize(new Dimension(220, 28));
        genomeFilter.setMinimumSize(new Dimension(180, 28));
        genomeFilter.setAlignmentX(0.0F);
        genomeFilter.addKeyListener(new KeyAdapter() {
            @Override
            public void keyReleased(KeyEvent e) {
                genomeEntryKeyReleased(e);
            }
        });
        filterPanel.add(genomeFilter, new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0,
                GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                new Insets(0, 0, 0, 0), 0, 0));

        contentPanel.add(filterPanel);

        //---- genomeList ----
        genomeList = new JList();
        genomeList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        genomeList.addListSelectionListener(e -> {
            GenomeListItem item = genomeList.getSelectedValue();
            if (item != null) {
                configureDownloadButtons(item);
            }
        });

        JScrollPane scrollPane1 = new JScrollPane();
        scrollPane1.setViewportView(genomeList);

        contentPanel.add(scrollPane1);

        //---- downloadSequenceCB ----
        downloadSequenceCB = new JCheckBox("Download sequence");
        downloadSequenceCB.setAlignmentX(1.0F);
        downloadSequenceCB.setToolTipText("Download the full sequence data for the genome.  Note that these files can be  large (human is about 1 gigabyte).");
        downloadSequenceCB.setMaximumSize(new Dimension(1000, 23));
        downloadSequenceCB.setPreferredSize(new Dimension(300, 23));
        downloadSequenceCB.setMinimumSize(new Dimension(300, 23));
        downloadSequenceCB.setEnabled(false);   // Disabled until a genome is selected.
        contentPanel.add(downloadSequenceCB);

        //---- downloadAnnotationsCB ----
        downloadAnnotationsCB = new JCheckBox("Download annotations");
        downloadAnnotationsCB.setAlignmentX(1.0F);
        downloadAnnotationsCB.setToolTipText("Download all annotation files referenced by the genome definition.");
        downloadAnnotationsCB.setMaximumSize(new Dimension(1000, 23));
        downloadAnnotationsCB.setPreferredSize(new Dimension(300, 23));
        downloadAnnotationsCB.setMinimumSize(new Dimension(300, 23));
        downloadAnnotationsCB.setEnabled(false);  // Disabled until a genome is selected
        contentPanel.add(downloadAnnotationsCB);


        dialogPane.add(contentPanel, BorderLayout.CENTER);

        //======== buttonBar ========
        JPanel buttonBar = new JPanel();
        buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
        buttonBar.setLayout(new GridBagLayout());
        ((GridBagLayout) buttonBar.getLayout()).columnWidths = new int[]{0, 85, 80};
        ((GridBagLayout) buttonBar.getLayout()).columnWeights = new double[]{1.0, 0.0, 0.0};

        //---- okButton ----
        JButton okButton = new JButton("OK");
        okButton.addActionListener(e -> okButtonActionPerformed(e));
        buttonBar.add(okButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
                GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                new Insets(0, 0, 5, 5), 0, 0));

        //---- cancelButton ---
        JButton cancelButton = new JButton("Cancel");
        cancelButton.addActionListener(e -> cancelButtonActionPerformed(e));
        buttonBar.add(cancelButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0,
                GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                new Insets(0, 0, 5, 0), 0, 0));

        dialogPane.add(buttonBar, BorderLayout.SOUTH);

        contentPane.add(dialogPane, BorderLayout.CENTER);
        pack();
        setLocationRelativeTo(getOwner());
    }


}
