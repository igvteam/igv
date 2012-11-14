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
import java.util.List;

/**
 * @author Jacob Silterra
 */
public class ManageGenomesDialog extends JDialog {

    private List<GenomeListItem> allListItems;
    private boolean cancelled = true;
    private List<GenomeListItem> removedValuesList = new ArrayList<GenomeListItem>();
    private GenomeListItem currentGenomeItem = null;

    public ManageGenomesDialog(Frame owner) {
        super(owner);
        initComponents();

        initData();

        genomeList.setCellRenderer(new ListCellRenderer() {
            @Override
            public Component getListCellRendererComponent(JList jList, Object value, int index, boolean isSelected, boolean cellHasFocus) {
                JLabel comp = new JLabel(value.toString());
                if (value instanceof GenomeListItem) {
                    GenomeListItem item = (GenomeListItem) value;
                    comp.setText(item.getDisplayableName());
                    comp.setToolTipText(item.getLocation());
                    if (isSelected) {
                        comp.setBackground(genomeList.getSelectionBackground());
                        comp.setForeground(genomeList.getSelectionForeground());
                        comp.setOpaque(isSelected);
                    }
                }
                return comp;
            }
        });
    }

    private void initData() {
        allListItems = new ArrayList<GenomeListItem>(GenomeManager.getInstance().getGenomes());
        String genomeId = GenomeManager.getInstance().getGenomeId();
        currentGenomeItem = GenomeManager.getInstance().getLoadedGenomeListItemById(genomeId);
        buildList();
        genomeList.setTransferHandler(new SimpleTransferHandler());
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
        cancelled = false;
        PreferenceManager.getInstance().saveGenomeIdDisplayList(allListItems);
        setVisible(false);
    }

    private void removeSelected() {
        List<GenomeListItem> selectedValuesList = genomeList.getSelectedValuesList();
        if (selectedValuesList.contains(currentGenomeItem)) {
            MessageUtils.showMessage("Cannot remove currently selected genome " + currentGenomeItem.getDisplayableName());
            return;
        }
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
        GenomeSelectionDialog dialog = new GenomeSelectionDialog(null, ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
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
        panel1 = new JPanel();
        addRemBar = new JPanel();
        addButton = new JButton();
        removeButton = new JButton();
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
    private JPanel panel1;
    private JPanel addRemBar;
    private JButton addButton;
    private JButton removeButton;
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

}
