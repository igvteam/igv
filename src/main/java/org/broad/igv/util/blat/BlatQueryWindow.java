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
 * Created by JFormDesigner on Fri Nov 30 11:55:45 EST 2012
 */

package org.broad.igv.util.blat;

import org.broad.igv.DirectoryManager;
import org.broad.igv.feature.PSLRecord;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.LongRunningTask;

import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.io.IOException;
import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

/**
 * @author Jim Robinson
 */
public class BlatQueryWindow extends JFrame {

    BlatTableModel model;

    public BlatQueryWindow(Component parent, String querySequence, java.util.List<PSLRecord> records) {

        if(parent != null) this.setLocationRelativeTo(parent);

        model = new BlatTableModel(records);
        initComponents();
        addSelectionListener();

        querySeqTextPane.setContentType("text/html");
        StringBuffer headerBuffer = new StringBuffer("<html>");
        headerBuffer.append("&nbsp;&nbsp;BLAT result for query sequence: <br>&nbsp;&nbsp&nbsp;&nbsp;");
        headerBuffer.append(querySequence);
        headerBuffer.append("<br><br>&nbsp;&nbsp;<i>Click on a row to go to alignment");
        querySeqTextPane.setText(headerBuffer.toString());
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
                    final String chr = model.getChr(selectedRow);
                    int start = model.getStart(selectedRow);
                    int end = model.getEnd(selectedRow);

                    // Expand region slightly for context
                    int w = (end - start) / 4;
                    final int estart = Math.max(0, start - w);
                    final int eend = end + w;

                    LongRunningTask.submit(new Runnable() {
                        @Override
                        public void run() {
                            IGV.getInstance().goToLocus(chr + ":" + estart + "-" + eend);
                        }
                    });
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
        headerPanel = new JPanel();
        panel1 = new JPanel();
        scrollPane2 = new JScrollPane();
        querySeqTextPane = new JTextPane();

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

            //======== headerPanel ========
            {
                headerPanel.setLayout(new BorderLayout());

                //======== panel1 ========
                {
                    panel1.setLayout(new FlowLayout(FlowLayout.LEFT));
                }
                headerPanel.add(panel1, BorderLayout.NORTH);

                //======== scrollPane2 ========
                {
                    scrollPane2.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
                    scrollPane2.setViewportView(querySeqTextPane);
                }
                headerPanel.add(scrollPane2, BorderLayout.CENTER);
            }
            contentPanel.add(headerPanel, BorderLayout.NORTH);
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
    private JPanel headerPanel;
    private JPanel panel1;
    private JScrollPane scrollPane2;
    private JTextPane querySeqTextPane;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
