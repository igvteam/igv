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

package org.broad.igv.ui.commandbar;

import org.broad.igv.event.GenomeResetEvent;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.feature.genome.DotGenomeUtils;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.genome.GenomeListItem;
import org.broad.igv.ui.genome.GenomeListManager;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.LongRunningTask;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;


public class RemoveGenomesDialog extends org.broad.igv.ui.IGVDialog  {

    public static final String LOCAL_SEQUENCE_CHAR = "\u002A";
    private static Logger log = LogManager.getLogger(RemoveGenomesDialog.class);

    private List<GenomeListItem> allListItems;

    private boolean haveLocalSequence = false;

    public RemoveGenomesDialog(Frame owner) {
        super(owner);
        initComponents();

        initData();

        genomeList.setCellRenderer(new GenomeCellRenderer());
    }

    private void initData() {

        allListItems = new ArrayList<>(GenomeListManager.getInstance().getGenomeTableRecords());

        for (GenomeListItem item : allListItems) {
            if (DotGenomeUtils.getLocalFasta(item.getId()) != null) {
                haveLocalSequence = true;
                break;
            }
        }

        buildList();

        label2.setVisible(haveLocalSequence);
    }

    private void buildList() {
        String currentId = GenomeManager.getInstance().getGenomeId();
        List<GenomeListItem> filteredList = allListItems.stream()
                .filter((item) -> !item.getId().equals(currentId))
                .collect(Collectors.toList());
        genomeList.setListData(filteredList.toArray(new GenomeListItem[0]));
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        setVisible(false);
    }

    private void saveButtonActionPerformed(ActionEvent event) {

        Runnable runnable = () -> {
            List<GenomeListItem> selectedValuesList = genomeList.getSelectedValuesList();
            if (selectedValuesList != null && !selectedValuesList.isEmpty()) {

                // Remove from the dropdown list
                GenomeListManager.getInstance().removeItems(selectedValuesList);

                // Remove downloaded genomes, if any, and associated files
                try {
                    GenomeManager.getInstance().deleteDownloadedGenomes(selectedValuesList);
                } catch (IOException e) {
                    log.error("Error deleting genome files", e);
                    MessageUtils.showErrorMessage("Error deleting genome files", e);
                }

                // If the last genome selected (DEFAULT_GENOME) was removed reset the key
                String lastGenomeKey = PreferencesManager.getPreferences().get(Constants.DEFAULT_GENOME);
                for (GenomeListItem item : selectedValuesList) {
                    if (lastGenomeKey.equals(item.getId())) {
                        PreferencesManager.getPreferences().remove(Constants.DEFAULT_GENOME);
                        break;
                    }
                }
            }
            if (selectedValuesList.size() > 0) {
                IGVEventBus.getInstance().post(new GenomeResetEvent());
            }

        };

        LongRunningTask.submit(runnable);
        setVisible(false);
    }


    private void removeSelected() {
        List<GenomeListItem> selectedValuesList = genomeList.getSelectedValuesList();
        allListItems.removeAll(selectedValuesList);
        buildList();
    }

    private void genomeListKeyReleased(KeyEvent e) {
        if (e.getKeyCode() == KeyEvent.VK_DELETE || e.getKeyCode() == KeyEvent.VK_BACK_SPACE) {
            removeSelected();
        }
    }

    private void removeButtonActionPerformed(ActionEvent e) {
        removeSelected();
    }


    private class GenomeCellRenderer implements ListCellRenderer {
        @Override
        public Component getListCellRendererComponent(JList list, Object value, int index, boolean isSelected, boolean cellHasFocus) {

            JLabel comp = new JLabel(value.toString());

            GenomeListItem item = (GenomeListItem) value;
            String displayableName = item.getDisplayableName();

            comp.setToolTipText(item.getPath());
            if (isSelected) {
                comp.setBackground(genomeList.getSelectionBackground());
                comp.setForeground(genomeList.getSelectionForeground());
                comp.setOpaque(isSelected);
            }

            if (DotGenomeUtils.getLocalFasta(item.getId()) != null) {
                displayableName += LOCAL_SEQUENCE_CHAR;
            }

            comp.setText(displayableName);
            return comp;
        }
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        label1 = new JTextArea();
        contentPanel = new JPanel();
        scrollPane1 = new JScrollPane();
        genomeList = new JList<>();
        label2 = new JLabel();
        panel1 = new JPanel();
        addRemBar = new JPanel();
        separator1 = new JSeparator();
        buttonBar = new JPanel();
        okButton = new JButton();
        cancelButton = new JButton();

        //======== this ========
        setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        setModalityType(Dialog.ModalityType.DOCUMENT_MODAL);
        setTitle("Remove Genomes");
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setPreferredSize(new Dimension(270, 400));
            dialogPane.setLayout(new BorderLayout());

            //---- label1 ----
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

                }
                panel1.add(addRemBar);
                panel1.add(separator1);

                //======== buttonBar ========
                {
                    buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                    buttonBar.setPreferredSize(new Dimension(196, 51));
                    buttonBar.setLayout(new FlowLayout(FlowLayout.TRAILING));

                    //---- okButton ----
                    okButton.setText("Remove");
                    okButton.setMaximumSize(new Dimension(93, 29));
                    okButton.setMinimumSize(new Dimension(93, 29));
                    okButton.setPreferredSize(new Dimension(93, 29));
                    okButton.addActionListener(e -> saveButtonActionPerformed(e));
                    buttonBar.add(okButton);

                    //---- cancelButton ----
                    cancelButton.setText("Cancel");
                    cancelButton.setMinimumSize(new Dimension(93, 29));
                    cancelButton.setPreferredSize(new Dimension(93, 29));
                    cancelButton.setMaximumSize(new Dimension(93, 29));
                    cancelButton.addActionListener(e -> cancelButtonActionPerformed(e));
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
    private JList<GenomeListItem> genomeList;
    private JLabel label2;
    private JPanel panel1;
    private JPanel addRemBar;
    private JSeparator separator1;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
