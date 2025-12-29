/*
 * Created by JFormDesigner on Fri Sep 07 13:22:05 EDT 2012
 */

package org.igv.ui.commandbar;

import org.igv.event.GenomeResetEvent;
import org.igv.event.IGVEventBus;
import org.igv.feature.genome.DotGenomeUtils;
import org.igv.feature.genome.GenomeManager;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.ui.genome.GenomeListItem;
import org.igv.ui.genome.GenomeListManager;
import org.igv.ui.util.MessageUtils;
import org.igv.util.LongRunningTask;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


public class RemoveGenomesDialog extends org.igv.ui.IGVDialog  {

    public static final String LOCAL_SEQUENCE_CHAR = "*";
    private static final Logger log = LogManager.getLogger(RemoveGenomesDialog.class);

    private List<GenomeListItem> allListItems;
    private boolean haveLocalSequence = false;

    private JList<GenomeListItem> genomeList;
    private JLabel label2;

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
                .toList();
        genomeList.setListData(filteredList.toArray(new GenomeListItem[0]));
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        setVisible(false);
    }

    private void saveButtonActionPerformed(ActionEvent event) {

        // Get selected values on EDT before submitting background task
        List<GenomeListItem> selectedValuesList = genomeList.getSelectedValuesList();

        Runnable runnable = () -> {
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
                if (lastGenomeKey != null) {
                    for (GenomeListItem item : selectedValuesList) {
                        if (lastGenomeKey.equals(item.getId())) {
                            PreferencesManager.getPreferences().remove(Constants.DEFAULT_GENOME);
                            break;
                        }
                    }
                }

                IGVEventBus.getInstance().post(new GenomeResetEvent());
            }
        };

        LongRunningTask.submit(runnable);
        setVisible(false);
    }





    private class GenomeCellRenderer implements ListCellRenderer<GenomeListItem> {
        @Override
        public Component getListCellRendererComponent(JList<? extends GenomeListItem> list, GenomeListItem value, int index, boolean isSelected, boolean cellHasFocus) {

            JLabel comp = new JLabel(value.toString());

            String displayableName = value.getDisplayableName();

            comp.setToolTipText(value.getPath());
            if (isSelected) {
                comp.setBackground(genomeList.getSelectionBackground());
                comp.setForeground(genomeList.getSelectionForeground());
                comp.setOpaque(isSelected);
            }

            if (DotGenomeUtils.getLocalFasta(value.getId()) != null) {
                displayableName += LOCAL_SEQUENCE_CHAR;
            }

            comp.setText(displayableName);
            return comp;
        }
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        JPanel dialogPane = new JPanel();
        JTextArea label1 = new JTextArea();
        JPanel contentPanel = new JPanel();
        JScrollPane scrollPane1 = new JScrollPane();
        genomeList = new JList<>();
        label2 = new JLabel();
        JPanel panel1 = new JPanel();
        JPanel addRemBar = new JPanel();
        JSeparator separator1 = new JSeparator();
        JPanel buttonBar = new JPanel();
        JButton okButton = new JButton();
        JButton cancelButton = new JButton();

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
                    okButton.addActionListener(this::saveButtonActionPerformed);
                    buttonBar.add(okButton);

                    //---- cancelButton ----
                    cancelButton.setText("Cancel");
                    cancelButton.setMinimumSize(new Dimension(93, 29));
                    cancelButton.setPreferredSize(new Dimension(93, 29));
                    cancelButton.setMaximumSize(new Dimension(93, 29));
                    cancelButton.addActionListener(this::cancelButtonActionPerformed);
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

}
