/*
 * Created by JFormDesigner on Fri Nov 30 11:55:45 EST 2012
 */

package org.broad.igv.blat;

import org.broad.igv.DirectoryManager;
import org.broad.igv.feature.PSLRecord;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.FileUtils;

import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.io.IOException;
import java.util.*;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

/**
 * @author Jim Robinson
 */
public class BlatResultsWindow extends JFrame {

    BlatTableModel model;

    public BlatResultsWindow(String querySequence, java.util.List<PSLRecord> records) {

        model = new BlatTableModel(records);
        initComponents();
        addSelectionListener();

        StringBuffer headerBuffer = new StringBuffer("<html>");
        headerBuffer.append("&nbsp;&nbsp;BLAT result for query sequence: <br>&nbsp;&nbsp&nbsp;&nbsp;");
        headerBuffer.append(querySequence);
        headerBuffer.append("<br><br>&nbsp;&nbsp;<i>Click on a row to go to alignment");
        headerLabel.setText(headerBuffer.toString());
    }

    private void addSelectionListener() {

        blatTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        ListSelectionModel rowSM = blatTable.getSelectionModel();
        rowSM.addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent e) {
                //Ignore extra messages.
                if (e.getValueIsAdjusting()) return;

                ListSelectionModel lsm = (ListSelectionModel) e.getSource();
                if (!lsm.isSelectionEmpty()) {
                    int selectedRow = lsm.getMinSelectionIndex();
                    String chr = model.getChr(selectedRow);
                    int start = model.getStart(selectedRow);
                    int end = model.getEnd(selectedRow);

                    // Expand region slightly for context
                    int w = (end - start) / 4;
                    start = Math.max(0, start - w);
                    end = end + w;

                    IGV.getInstance().goToLocus(chr + ":" + start + "-" + end);
                }
            }
        });
    }


    private void closeItemActionPerformed(ActionEvent e) {
        setVisible(false);
        dispose();
    }

    private void saveItemActionPerformed(ActionEvent e) {

        File f = FileDialogUtils.chooseFile("Save BLAT results", DirectoryManager.getUserDirectory(), FileDialogUtils.SAVE);
        if (f != null) {
            try {
                model.save(f);
            } catch (IOException e1) {
                MessageUtils.showErrorMessage("Error saving blat results ", e1);
            }
        }
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        menuBar = new JMenuBar();
        menu1 = new JMenu();
        saveItem = new JMenuItem();
        closeItem = new JMenuItem();
        contentPanel = new JPanel();
        scrollPane1 = new JScrollPane();
        blatTable = new JTable(model);
        headerLabel = new JLabel();

        //======== this ========
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== menuBar ========
        {

            //======== menu1 ========
            {
                menu1.setText("File");

                //---- saveItem ----
                saveItem.setText("Save...");
                saveItem.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        saveItemActionPerformed(e);
                    }
                });
                menu1.add(saveItem);

                //---- closeItem ----
                closeItem.setText("Close");
                closeItem.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        closeItemActionPerformed(e);
                    }
                });
                menu1.add(closeItem);
            }
            menuBar.add(menu1);
        }
        setJMenuBar(menuBar);

        //======== contentPanel ========
        {
            contentPanel.setLayout(new BorderLayout(0, 5));

            //======== scrollPane1 ========
            {
                scrollPane1.setViewportView(blatTable);
            }
            contentPanel.add(scrollPane1, BorderLayout.CENTER);
            contentPanel.add(headerLabel, BorderLayout.NORTH);
        }
        contentPane.add(contentPanel, BorderLayout.CENTER);
        setSize(870, 570);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JMenuBar menuBar;
    private JMenu menu1;
    private JMenuItem saveItem;
    private JMenuItem closeItem;
    private JPanel contentPanel;
    private JScrollPane scrollPane1;
    private JTable blatTable;
    private JLabel headerLabel;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
