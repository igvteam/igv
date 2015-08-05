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

package org.broad.igv.ui;

import org.broad.igv.DirectoryManager;
import org.broad.igv.feature.genome.GenomeListItem;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.RunnableResult;

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
public class GenomeSelectionDialog extends javax.swing.JDialog {

    private boolean isCanceled = true;
    private List<GenomeListItem> selectedValuesList = null;
    private List<GenomeListItem> allListItems;
    private DefaultListModel genomeListModel;

    /**
     * @param parent
     * @param listSelectionMode Selection mode for genome list
     */
    public GenomeSelectionDialog(java.awt.Frame parent, Collection<GenomeListItem> inputListItems, final int listSelectionMode) {
        super(parent);
        initComponents();
        UIUtilities.invokeOnEventThread(new Runnable() {

            @Override
            public void run() {
                genomeList.setSelectionMode(listSelectionMode);
            }
        });

        initData(inputListItems);

        downloadSequenceCB.setVisible(listSelectionMode == ListSelectionModel.SINGLE_SELECTION);
    }

    private void initData(Collection<GenomeListItem> inputListItems) {
        this.allListItems = new ArrayList<GenomeListItem>(inputListItems);
        //We don't show those currently in the box
        //We only check by id, as some local genomes might be the same except for path
        for(GenomeListItem gli: GenomeManager.getInstance().getGenomes()){
            removeById(gli.getId());
        }
        this.allListItems.removeAll(GenomeManager.getInstance().getGenomes());
        rebuildGenomeList();
    }

    private void removeById(String id){
        int ci = 0;
        for(; ci < this.allListItems.size(); ci++){
            GenomeListItem gli = this.allListItems.get(ci);
            if(gli.getId().equals(id)){
                break;
            }
        }
        if(ci < this.allListItems.size()) this.allListItems.remove(ci);
    }

    private void rebuildGenomeList() {
        String filterText = genomeFilter.getText().trim().toLowerCase();
        rebuildGenomeList(filterText);
    }

    /**
     * Filter the list of displayed genomes so we only show this
     * with the text the user entered.
     */
    private void rebuildGenomeList(String filterText) {
        if (genomeListModel == null) {
            genomeListModel = new DefaultListModel();
            UIUtilities.invokeOnEventThread(new Runnable() {
                @Override
                public void run() {
                    genomeList.setModel(genomeListModel);
                }
            });
        }
        genomeListModel.clear();
        filterText = filterText.toLowerCase().trim();
        for (GenomeListItem listItem : allListItems) {
            if (listItem.getDisplayableName().toLowerCase().contains(filterText)) {
                genomeListModel.addElement(listItem);
            }
        }
    }

    /**
     * If a genome is single clicked, we just store the selection.
     * When a genome is double clicked, we treat that as the user
     * wanting to load the genome.
     *
     * @param e
     */
    private void genomeListMouseClicked(MouseEvent e) {
        switch (e.getClickCount()) {
            case 1:
                List<GenomeListItem> selValues = genomeList.getSelectedValuesList();
                downloadSequenceCB.setEnabled(selValues != null && selValues.size() == 1);
                break;
            case 2:
                okButtonActionPerformed(null);
                break;
        }
    }

    private void genomeEntryKeyReleased(KeyEvent e) {
        rebuildGenomeList(genomeFilter.getText());
    }

    public List<GenomeListItem> getSelectedValuesList() {
        return selectedValuesList;
    }

    public boolean downloadSequence(){
        return !isCanceled() && downloadSequenceCB.isEnabled() && downloadSequenceCB.isSelected();
    }

    static void downloadGenome(final Frame dialogsParent, final GenomeListItem genomeListItem) {
        final File targetDir = DirectoryManager.getGenomeCacheDirectory();
        if(!targetDir.exists()){
            targetDir.mkdirs();
        }

        Runnable runnable = new Runnable() {
            @Override
            public void run() {

                try {
                    RunnableResult result = GenomeManager.getInstance().downloadWholeGenome(genomeListItem.getLocation(), targetDir, dialogsParent);
                    if(result == RunnableResult.FAILURE){
                        throw new IOException("Unknown Failure");
                    }
                } catch (IOException e) {
                    String msg = String.format("Error downloading genome %s from %s: %s", genomeListItem.getId(), genomeListItem.getLocation(), e.getMessage());
                    MessageUtils.showErrorMessage(msg, e);
                }
            }
        };

        LongRunningTask.submit(runnable);
    }

    public boolean isCanceled() {
        return isCanceled;
    }

    /**
     * This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    // Generated using JFormDesigner non-commercial license
    private void initComponents() {
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        textArea1 = new JTextArea();
        filterPanel = new JPanel();
        label1 = new JLabel();
        genomeFilter = new JTextField();
        scrollPane1 = new JScrollPane();
        genomeList = new JList7<GenomeListItem>();
        downloadSequenceCB = new JCheckBox();
        buttonBar = new JPanel();
        okButton = new JButton();
        cancelButton = new JButton();

        //======== this ========
        setModal(true);
        setTitle("Genomes to add to list");
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setPreferredSize(new Dimension(350, 500));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.Y_AXIS));

                //---- textArea1 ----
                textArea1.setText("Selected genomes will be added to the genome dropdown list.");
                textArea1.setLineWrap(true);
                textArea1.setWrapStyleWord(true);
                textArea1.setBackground(UIManager.getColor("Button.background"));
                textArea1.setRows(2);
                textArea1.setMaximumSize(new Dimension(2147483647, 60));
                textArea1.setRequestFocusEnabled(false);
                textArea1.setEditable(false);
                contentPanel.add(textArea1);

                //======== filterPanel ========
                {
                    filterPanel.setMaximumSize(new Dimension(2147483647, 28));
                    filterPanel.setLayout(new GridBagLayout());
                    ((GridBagLayout)filterPanel.getLayout()).columnWidths = new int[] {0, 0, 0};
                    ((GridBagLayout)filterPanel.getLayout()).rowHeights = new int[] {0, 0};
                    ((GridBagLayout)filterPanel.getLayout()).columnWeights = new double[] {1.0, 1.0, 1.0E-4};
                    ((GridBagLayout)filterPanel.getLayout()).rowWeights = new double[] {1.0, 1.0E-4};

                    //---- label1 ----
                    label1.setText("Filter:");
                    label1.setLabelFor(genomeFilter);
                    label1.setRequestFocusEnabled(false);
                    filterPanel.add(label1, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0,
                        GridBagConstraints.WEST, GridBagConstraints.VERTICAL,
                        new Insets(0, 0, 0, 0), 0, 0));

                    //---- genomeFilter ----
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
                }
                contentPanel.add(filterPanel);

                //======== scrollPane1 ========
                {

                    //---- genomeList ----
                    genomeList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
                    genomeList.addMouseListener(new MouseAdapter() {
                        @Override
                        public void mouseClicked(MouseEvent e) {
                            genomeListMouseClicked(e);
                        }
                    });
                    scrollPane1.setViewportView(genomeList);
                }
                contentPanel.add(scrollPane1);

                //---- downloadSequenceCB ----
                downloadSequenceCB.setText("Download Sequence");
                downloadSequenceCB.setAlignmentX(1.0F);
                downloadSequenceCB.setToolTipText("Download the full sequence for this organism. Note that these files can be very large (human is about 3 gigabytes)");
                downloadSequenceCB.setMaximumSize(new Dimension(1000, 23));
                downloadSequenceCB.setPreferredSize(new Dimension(300, 23));
                downloadSequenceCB.setMinimumSize(new Dimension(300, 23));
                contentPanel.add(downloadSequenceCB);
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
                okButton.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        okButtonActionPerformed(e);
                    }
                });
                buttonBar.add(okButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 5, 5), 0, 0));

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
                    new Insets(0, 0, 5, 0), 0, 0));
            }
            dialogPane.add(buttonBar, BorderLayout.SOUTH);
        }
        contentPane.add(dialogPane, BorderLayout.CENTER);
        pack();
        setLocationRelativeTo(getOwner());
    }// </editor-fold>//GEN-END:initComponents

    private void cancelButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_cancelButtonActionPerformed
        isCanceled = true;
        selectedValuesList = null;
        setVisible(false);
        dispose();
    }//GEN-LAST:event_cancelButtonActionPerformed

    private void okButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_okButtonActionPerformed
        isCanceled = false;
        selectedValuesList = genomeList.getSelectedValuesList();
        setVisible(false);
        dispose();
    }//GEN-LAST:event_okButtonActionPerformed

    // Variables declaration - do not modify//GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JTextArea textArea1;
    private JPanel filterPanel;
    private JLabel label1;
    private JTextField genomeFilter;
    private JScrollPane scrollPane1;
    private JList7<GenomeListItem> genomeList;
    private JCheckBox downloadSequenceCB;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    // End of variables declaration//GEN-END:variables

}
