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
 * Created by JFormDesigner on Fri Sep 07 13:22:05 EDT 2012
 */

package org.broad.igv.ui;

import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.genome.GenomeListItem;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.util.MessageUtils;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * @author Jacob Silterra
 */
public class ManageGenomesDialog extends JDialog {

    private List<GenomeListItem> allListItems;
    private boolean cancelled = true;
    private List<GenomeListItem> removedValuesList = new ArrayList<GenomeListItem>();
    private GenomeListItem currentGenomeItem = null;

    private boolean haveLocalGenomes = false;
    private static final String LOCAL_SEQUENCE_CHAR = "\u002A";

    public ManageGenomesDialog(Frame owner) {
        super(owner);
        initComponents();

        initData();

        genomeList.setCellRenderer(new GenomeCellRenderer());
    }

    private void initData() {
        allListItems = new ArrayList<GenomeListItem>(GenomeManager.getInstance().getGenomes());
        for(GenomeListItem item: allListItems){
            if(item.hasDownloadedSequence()){
                haveLocalGenomes = true;
                break;
            }
        }
        String genomeId = GenomeManager.getInstance().getGenomeId();
        currentGenomeItem = GenomeManager.getInstance().getLoadedGenomeListItemById(genomeId);
        buildList();
        genomeList.setTransferHandler(new SimpleTransferHandler());

        addButton.setEnabled(!GenomeManager.getInstance().isServerGenomeListUnreachable());
        label2.setVisible(haveLocalGenomes);
    }


    private void buildList() {
        genomeList.setListData(allListItems.toArray());
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        cancelled = true;
        removedValuesList = null;
        setVisible(false);
    }

    private void saveButtonActionPerformed(ActionEvent e) {
        if(removedValuesList.size() > 0){
            int countHasSeq = 0;
            for(GenomeListItem item: removedValuesList){
                if(item.hasDownloadedSequence()){
                    countHasSeq++;
                }
            }

            if(countHasSeq > 0){
                String msg = String.format("%d of the genomes you chose to remove have downloaded sequences. Those will be deleted as well. Are you sure?", countHasSeq);
                boolean sure = MessageUtils.confirm(msg);
                if(!sure) return;
            }
        }

        cancelled = false;
        PreferenceManager.getInstance().saveGenomeIdDisplayList(allListItems);
        setVisible(false);
    }

    private void removeSelected() {
        List<GenomeListItem> selectedValuesList = genomeList.getSelectedValuesList();
//        if (selectedValuesList.contains(currentGenomeItem)) {
//            MessageUtils.showMessage("Cannot remove currently selected genome " + currentGenomeItem.getDisplayableName());
//            return;
//        }
        removedValuesList.addAll(selectedValuesList);
        allListItems.removeAll(selectedValuesList);
        buildList();
    }

    public List<GenomeListItem> getRemovedValuesList() {
        return removedValuesList;
    }

    private void genomeListKeyReleased(KeyEvent e) {
        if (e.getKeyCode() == KeyEvent.VK_DELETE || e.getKeyCode() == KeyEvent.VK_BACK_SPACE) {
            removeSelected();
        }
    }

    private void removeButtonActionPerformed(ActionEvent e) {
        removeSelected();
    }

    private void addButtonActionPerformed(ActionEvent e) {
        Collection<GenomeListItem> inputListItems = GenomeManager.getInstance().getGenomeArchiveList();
        if(inputListItems == null){
            IOException exc = new IOException("Unable to reach genome server");
            MessageUtils.showErrorMessage(exc.getMessage(), exc);
            return;
        }
        GenomeSelectionDialog dialog = new GenomeSelectionDialog(null, inputListItems, ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        dialog.setVisible(true);
        List<GenomeListItem> selectedValues = dialog.getSelectedValuesList();
        if (selectedValues != null) {
            allListItems.addAll(selectedValues);
            buildList();
        }
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        label1 = new JTextArea();
        contentPanel = new JPanel();
        scrollPane1 = new JScrollPane();
        genomeList = new JList7<GenomeListItem>();
        label2 = new JLabel();
        panel1 = new JPanel();
        addRemBar = new JPanel();
        addButton = new JButton();
        removeButton = new JButton();
        separator1 = new JSeparator();
        buttonBar = new JPanel();
        okButton = new JButton();
        cancelButton = new JButton();

        //======== this ========
        setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        setModalityType(Dialog.ModalityType.DOCUMENT_MODAL);
        setTitle("Manage Genome List");
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setPreferredSize(new Dimension(270, 400));
            dialogPane.setLayout(new BorderLayout());

            //---- label1 ----
            label1.setText("Drag and drop genomes to change their order in the genome list. \nSelect and press delete, or click \"Remove\", to remove them.");
            label1.setRows(2);
            label1.setEditable(false);
            label1.setBackground(UIManager.getColor("Button.background"));
            label1.setWrapStyleWord(true);
            label1.setLineWrap(true);
            dialogPane.add(label1, BorderLayout.NORTH);

            //======== contentPanel ========
            {
                contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.Y_AXIS));

                //======== scrollPane1 ========
                {

                    //---- genomeList ----
                    genomeList.setMaximumSize(new Dimension(39, 5000));
                    genomeList.setDropMode(DropMode.INSERT);
                    genomeList.setDragEnabled(true);
                    genomeList.addKeyListener(new KeyAdapter() {
                        @Override
                        public void keyReleased(KeyEvent e) {
                            genomeListKeyReleased(e);
                        }
                    });
                    scrollPane1.setViewportView(genomeList);
                }
                contentPanel.add(scrollPane1);

                //---- label2 ----
                label2.setText("Sequence on local machine");
                label2.setLabelFor(genomeList);
                label2.setAlignmentX(1.0F);
                label2.setComponentOrientation(ComponentOrientation.LEFT_TO_RIGHT);
                label2.setPreferredSize(new Dimension(400, 16));
                label2.setMaximumSize(new Dimension(400, 16));
                label2.setMinimumSize(new Dimension(100, 16));
                label2.setText(LOCAL_SEQUENCE_CHAR + label2.getText());
                contentPanel.add(label2);
            }
            dialogPane.add(contentPanel, BorderLayout.CENTER);

            //======== panel1 ========
            {
                panel1.setLayout(new BoxLayout(panel1, BoxLayout.Y_AXIS));

                //======== addRemBar ========
                {
                    addRemBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                    addRemBar.setPreferredSize(new Dimension(196, 51));
                    addRemBar.setMinimumSize(new Dimension(201, 51));
                    addRemBar.setLayout(new FlowLayout(FlowLayout.TRAILING, 1, 5));

                    //---- addButton ----
                    addButton.setText("Add From Server");
                    addButton.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            addButtonActionPerformed(e);
                        }
                    });
                    addRemBar.add(addButton);

                    //---- removeButton ----
                    removeButton.setText("Remove");
                    removeButton.setToolTipText("Remove selected genomes from list");
                    removeButton.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            removeButtonActionPerformed(e);
                        }
                    });
                    addRemBar.add(removeButton);
                }
                panel1.add(addRemBar);
                panel1.add(separator1);

                //======== buttonBar ========
                {
                    buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                    buttonBar.setPreferredSize(new Dimension(196, 51));
                    buttonBar.setLayout(new FlowLayout(FlowLayout.TRAILING));

                    //---- okButton ----
                    okButton.setText("Save");
                    okButton.setMaximumSize(new Dimension(93, 29));
                    okButton.setMinimumSize(new Dimension(93, 29));
                    okButton.setPreferredSize(new Dimension(93, 29));
                    okButton.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            saveButtonActionPerformed(e);
                        }
                    });
                    buttonBar.add(okButton);

                    //---- cancelButton ----
                    cancelButton.setText("Cancel");
                    cancelButton.setMinimumSize(new Dimension(93, 29));
                    cancelButton.setPreferredSize(new Dimension(93, 29));
                    cancelButton.setMaximumSize(new Dimension(93, 29));
                    cancelButton.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            cancelButtonActionPerformed(e);
                        }
                    });
                    buttonBar.add(cancelButton);
                }
                panel1.add(buttonBar);
            }
            dialogPane.add(panel1, BorderLayout.SOUTH);
        }
        contentPane.add(dialogPane, BorderLayout.CENTER);
        pack();
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JTextArea label1;
    private JPanel contentPanel;
    private JScrollPane scrollPane1;
    private JList7<GenomeListItem> genomeList;
    private JLabel label2;
    private JPanel panel1;
    private JPanel addRemBar;
    private JButton addButton;
    private JButton removeButton;
    private JSeparator separator1;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables

    public boolean isCancelled() {
        return cancelled;
    }

    private int findItem(String text) {
        int index = 0;
        for (GenomeListItem item : allListItems) {
            if (item.getId().equals(text) || item.getDisplayableName().equals(text)) {
                return index;
            }
            index++;
        }
        return -1;

    }

    private class SimpleTransferHandler extends TransferHandler {
        @Override
        public int getSourceActions(JComponent c) {
            return TransferHandler.MOVE;
        }

        @Override
        protected Transferable createTransferable(JComponent c) {
            return new StringSelection(PreferenceManager.generateGenomeIdString(genomeList.getSelectedValuesList()));
        }

        @Override
        public boolean importData(TransferSupport support) {
            if (!canImport(support)) {
                return false;
            }
            JList.DropLocation dropLocation = (JList.DropLocation) support.getDropLocation();
            int toIndex = dropLocation.getIndex();
            String[] genomeIds;
            try {
                String genomeIdString = (String) support.getTransferable().getTransferData(DataFlavor.stringFlavor);
                genomeIds = genomeIdString.split(PreferenceManager.HISTORY_DELIMITER);
            } catch (UnsupportedFlavorException e) {
                return false;
            } catch (IOException e) {
                return false;
            }

            if (genomeIds == null || genomeIds.length == 0) {
                return false;
            }

            int numMoved = 0;
            for (String genomeId : genomeIds) {
                int fromIndex = findItem(genomeId);
                if (fromIndex < 0 || fromIndex >= allListItems.size() || fromIndex == toIndex) {
                    continue;
                }
                //We need to account for the fact that the proper
                //insertion location is one smaller, once the item being moved
                //is removed.
                if (toIndex > fromIndex) toIndex--;
                GenomeListItem item = allListItems.remove(fromIndex);
                allListItems.add(toIndex, item);
                numMoved++;
                //Account for adding multiple items, want to add them to successive indices
                toIndex++;
            }
            buildList();
            return numMoved > 0;
        }

        @Override
        public boolean canImport(TransferSupport support) {
            support.setShowDropLocation(true);
            return support.isDrop();
        }
    }

    private class GenomeCellRenderer implements ListCellRenderer{
        @Override
        public Component getListCellRendererComponent(JList list, Object value, int index, boolean isSelected, boolean cellHasFocus) {

            JLabel comp = new JLabel(value.toString());

            GenomeListItem item = (GenomeListItem) value;
            String displayableName = item.getDisplayableName();

            comp.setToolTipText(item.getLocation());
            if (isSelected) {
                comp.setBackground(genomeList.getSelectionBackground());
                comp.setForeground(genomeList.getSelectionForeground());
                comp.setOpaque(isSelected);
            }

            if(item.hasDownloadedSequence()){
                displayableName += LOCAL_SEQUENCE_CHAR;
            }

            comp.setText(displayableName);
            return comp;
        }
    }

}
