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
 * @author User #2
 */
public class ManageGenomesDialog extends JDialog {

    private List<GenomeListItem> allListItems;
    private boolean cancelled = true;

    public ManageGenomesDialog(Frame owner) {
        super(owner);
        initComponents();

        initData();
    }

    private void initData() {
        allListItems = new ArrayList<GenomeListItem>(GenomeManager.getInstance().getGenomes());
        buildList();
        genomeList.setTransferHandler(new SimpleTransferHandler());
    }


    private void buildList() {
        genomeList.setListData(allListItems.toArray());
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        cancelled = true;
        setVisible(false);
    }

    private void okButtonActionPerformed(ActionEvent e) {
        cancelled = false;
        PreferenceManager.getInstance().saveGenomeIdDisplayList(allListItems);
        setVisible(false);
    }

    private void genomeListKeyReleased(KeyEvent e) {
        if (e.getKeyCode() == KeyEvent.VK_DELETE || e.getKeyCode() == KeyEvent.VK_BACK_SPACE) {
            allListItems.removeAll(genomeList.getSelectedValuesList());
            buildList();
        }
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        scrollPane1 = new JScrollPane();
        genomeList = new JList();
        buttonBar = new JPanel();
        okButton = new JButton();
        cancelButton = new JButton();

        //======== this ========
        setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        setModalityType(Dialog.ModalityType.DOCUMENT_MODAL);
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setPreferredSize(new Dimension(270, 400));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.Y_AXIS));

                //======== scrollPane1 ========
                {

                    //---- genomeList ----
                    genomeList.setMaximumSize(new Dimension(39, 5000));
                    genomeList.setDropMode(DropMode.INSERT);
                    genomeList.setDragEnabled(true);
                    genomeList.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);
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

            //======== buttonBar ========
            {
                buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                buttonBar.setPreferredSize(new Dimension(196, 51));
                buttonBar.setLayout(new GridBagLayout());
                ((GridBagLayout) buttonBar.getLayout()).columnWidths = new int[]{0, 85, 80};
                ((GridBagLayout) buttonBar.getLayout()).columnWeights = new double[]{1.0, 0.0, 0.0};

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
        pack();
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JScrollPane scrollPane1;
    private JList genomeList;
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
