/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

/*
 * Created by JFormDesigner on Mon Dec 13 14:39:17 EST 2010
 */

package org.broad.igv.lists;

import java.awt.*;
import java.awt.event.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Map;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;

import com.jidesoft.swing.*;

/**
 * @author Stan Diamond
 */
public class GeneListManagerUI extends JDialog {

    ListListModel listModel;
    GeneListModel geneListModel;
    Map<String, GeneList> geneLists;

    public GeneListManagerUI(Frame owner) {
        super(owner);
        initComponents();

        geneLists = GeneListManager.getGeneLists();

        groups.setModel(new AbstractListModel() {
            String[] values = {
                    "All"
            };

            public int getSize() {
                return values.length;
            }

            public Object getElementAt(int i) {
                return values[i];
            }
        });

        listModel = new ListListModel();
        lists.setModel(listModel);

        geneListModel = new GeneListModel();
        loci.setModel(geneListModel);
    }


    private void groupsValueChanged(ListSelectionEvent e) {
        System.out.println(e.getFirstIndex());
    }

    private void listsValueChanged(ListSelectionEvent e) {
        String listName = (String) lists.getSelectedValue();
        if (listName == null) {
            geneListModel.clear();
        } else {
            GeneList gl = listModel.getGeneList(listName);
            geneListModel.setGeneList(gl);
            loci.setModel(geneListModel);
            loci.updateUI();
        }
    }

    private void listLabelMouseClicked(MouseEvent e) {
        listModel.sort();
        lists.updateUI();
    }

    private void searchBoxKeyReleased(KeyEvent e) {

        listModel.filter();
        lists.clearSelection();
        lists.updateUI();
        loci.updateUI();
    }


    private void editButtonActionPerformed(ActionEvent e) {

        String selection = (String) lists.getSelectedValue();
        if (selection != null) {
            GeneList geneList = geneLists.get(selection);
            GeneListInputDialog dlg = new GeneListInputDialog(this, geneList);
            dlg.setVisible(true);
        }
    }


    class ListListModel extends AbstractListModel {

        boolean sortAscending;
        ArrayList<String> filteredNames;

        ListListModel() {
            filteredNames = new ArrayList(geneLists.size());
            filter();
        }

        public int getSize() {
            return filteredNames.size();
        }

        public Object getElementAt(int i) {
            return filteredNames.get(i);
        }

        public int getIndexFor(String str) {
            return filteredNames.indexOf(str);
        }

        void sort() {

            Collections.sort(filteredNames, new Comparator<String>() {
                public int compare(String s, String s1) {
                    return sortAscending ? s.compareTo(s1) : s1.compareTo(s);
                }
            });
            sortAscending = !sortAscending;
        }

        GeneList getGeneList(String listName) {
            return geneLists.get(listName);
        }

        void filter() {
            filteredNames.clear();
            for (Map.Entry<String, GeneList> entry : geneLists.entrySet()) {
                String name = entry.getKey();
                GeneList gl = entry.getValue();

                if (gl != null) {
                    if (gl != null && isPassFilter(gl)) {
                        filteredNames.add(name);
                    }
                }
            }
        }

        boolean isPassFilter(GeneList geneList) {
            String filterString = searchBox.getText();
            if (filterString != null && filterString.trim().length() > 0) {
                String tmp = filterString.trim().toLowerCase();

                if (geneList.getName().toLowerCase().contains(tmp)) {
                    return true;
                }

                for (String gene : geneList.getLoci()) {
                    if (gene.toLowerCase().contains(tmp)) {
                        return true;
                    }
                }

                return false;

            }
            return true;
        }
    }


    class GeneListModel extends AbstractListModel {

        java.util.List<String> genes;

        GeneListModel() {
            genes = new ArrayList();
        }

        void setGeneList(GeneList geneList) {

            genes = geneList == null ? new ArrayList() : geneList.getLoci();
        }

        public int getSize() {
            return genes.size();
        }

        public Object getElementAt(int i) {
            return genes.get(i);
        }

        public void clear() {
            genes.clear();
        }
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        panel2 = new JPanel();
        panel1 = new JPanel();
        panel6 = new JPanel();
        label1 = new JLabel();
        searchBox = new JTextField();
        splitPane2 = new JSplitPane();
        splitPane1 = new JSplitPane();
        panel3 = new JPanel();
        label3 = new JLabel();
        scrollPane1 = new JScrollPane();
        groups = new JList();
        panel7 = new JPanel();
        addGroupButton = new JideButton();
        panel4 = new JPanel();
        listLabel = new JLabel();
        scrollPane2 = new JScrollPane();
        lists = new JList();
        panel8 = new JPanel();
        jideButton2 = new JideButton();
        editButton = new JideButton();
        panel5 = new JPanel();
        label4 = new JLabel();
        scrollPane3 = new JScrollPane();
        loci = new JList();
        panel9 = new JPanel();
        jideButton4 = new JideButton();
        buttonBar = new JPanel();
        okButton = new JButton();

        //======== this ========
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setBorder(new EtchedBorder(EtchedBorder.RAISED));
                contentPanel.setLayout(new BorderLayout());

                //======== panel2 ========
                {
                    panel2.setBackground(new Color(204, 204, 204));
                    panel2.setLayout(new BoxLayout(panel2, BoxLayout.X_AXIS));

                    //======== panel1 ========
                    {
                        panel1.setBorder(new EmptyBorder(5, 5, 5, 0));
                        panel1.setBackground(new Color(204, 204, 204));
                        panel1.setLayout(new BoxLayout(panel1, BoxLayout.X_AXIS));

                        //======== panel6 ========
                        {
                            panel6.setMinimumSize(new Dimension(400, 0));
                            panel6.setPreferredSize(new Dimension(400, 0));
                            panel6.setOpaque(false);
                            panel6.setLayout(null);
                        }
                        panel1.add(panel6);

                        //---- label1 ----
                        label1.setText("Search");
                        label1.setBorder(new EmptyBorder(0, 0, 0, 10));
                        panel1.add(label1);

                        //---- searchBox ----
                        searchBox.setBorder(new BevelBorder(BevelBorder.LOWERED));
                        searchBox.addKeyListener(new KeyAdapter() {
                            @Override
                            public void keyReleased(KeyEvent e) {
                                searchBoxKeyReleased(e);
                            }
                        });
                        panel1.add(searchBox);
                    }
                    panel2.add(panel1);
                }
                contentPanel.add(panel2, BorderLayout.NORTH);

                //======== splitPane2 ========
                {
                    splitPane2.setBorder(null);
                    splitPane2.setDividerLocation(400);

                    //======== splitPane1 ========
                    {
                        splitPane1.setBorder(null);
                        splitPane1.setDividerLocation(150);

                        //======== panel3 ========
                        {
                            panel3.setBorder(LineBorder.createBlackLineBorder());
                            panel3.setLayout(new BorderLayout());

                            //---- label3 ----
                            label3.setText("Group");
                            panel3.add(label3, BorderLayout.NORTH);

                            //======== scrollPane1 ========
                            {

                                //---- groups ----
                                groups.setModel(new AbstractListModel() {
                                    String[] values = {
                                        "All"
                                    };
                                    public int getSize() { return values.length; }
                                    public Object getElementAt(int i) { return values[i]; }
                                });
                                groups.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
                                groups.addListSelectionListener(new ListSelectionListener() {
                                    public void valueChanged(ListSelectionEvent e) {
                                        groupsValueChanged(e);
                                    }
                                });
                                scrollPane1.setViewportView(groups);
                            }
                            panel3.add(scrollPane1, BorderLayout.CENTER);

                            //======== panel7 ========
                            {
                                panel7.setLayout(new BoxLayout(panel7, BoxLayout.X_AXIS));

                                //---- addGroupButton ----
                                addGroupButton.setText("+");
                                panel7.add(addGroupButton);
                            }
                            panel3.add(panel7, BorderLayout.SOUTH);
                        }
                        splitPane1.setLeftComponent(panel3);

                        //======== panel4 ========
                        {
                            panel4.setBorder(LineBorder.createBlackLineBorder());
                            panel4.setLayout(new BorderLayout());

                            //---- listLabel ----
                            listLabel.setText("List");
                            listLabel.addMouseListener(new MouseAdapter() {
                                @Override
                                public void mouseClicked(MouseEvent e) {
                                    listLabelMouseClicked(e);
                                }
                            });
                            panel4.add(listLabel, BorderLayout.PAGE_START);

                            //======== scrollPane2 ========
                            {

                                //---- lists ----
                                lists.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
                                lists.addListSelectionListener(new ListSelectionListener() {
                                    public void valueChanged(ListSelectionEvent e) {
                                        listsValueChanged(e);
                                    }
                                });
                                scrollPane2.setViewportView(lists);
                            }
                            panel4.add(scrollPane2, BorderLayout.CENTER);

                            //======== panel8 ========
                            {
                                panel8.setLayout(new BoxLayout(panel8, BoxLayout.X_AXIS));

                                //---- jideButton2 ----
                                jideButton2.setText("+");
                                panel8.add(jideButton2);

                                //---- editButton ----
                                editButton.setText("Edit");
                                editButton.addActionListener(new ActionListener() {
                                    public void actionPerformed(ActionEvent e) {
                                        editButtonActionPerformed(e);
                                    }
                                });
                                panel8.add(editButton);
                            }
                            panel4.add(panel8, BorderLayout.PAGE_END);
                        }
                        splitPane1.setRightComponent(panel4);
                    }
                    splitPane2.setLeftComponent(splitPane1);

                    //======== panel5 ========
                    {
                        panel5.setBorder(LineBorder.createBlackLineBorder());
                        panel5.setLayout(new BorderLayout());

                        //---- label4 ----
                        label4.setText("Loci");
                        panel5.add(label4, BorderLayout.NORTH);

                        //======== scrollPane3 ========
                        {
                            scrollPane3.setViewportView(loci);
                        }
                        panel5.add(scrollPane3, BorderLayout.CENTER);

                        //======== panel9 ========
                        {
                            panel9.setLayout(new BoxLayout(panel9, BoxLayout.X_AXIS));

                            //---- jideButton4 ----
                            jideButton4.setText(" ");
                            panel9.add(jideButton4);
                        }
                        panel5.add(panel9, BorderLayout.SOUTH);
                    }
                    splitPane2.setRightComponent(panel5);
                }
                contentPanel.add(splitPane2, BorderLayout.CENTER);
            }
            dialogPane.add(contentPanel, BorderLayout.CENTER);

            //======== buttonBar ========
            {
                buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                buttonBar.setLayout(new GridBagLayout());
                ((GridBagLayout)buttonBar.getLayout()).columnWidths = new int[] {0, 80};
                ((GridBagLayout)buttonBar.getLayout()).columnWeights = new double[] {1.0, 0.0};

                //---- okButton ----
                okButton.setText("OK");
                buttonBar.add(okButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 0, 0), 0, 0));
            }
            dialogPane.add(buttonBar, BorderLayout.SOUTH);
        }
        contentPane.add(dialogPane, BorderLayout.CENTER);
        setSize(715, 575);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JPanel panel2;
    private JPanel panel1;
    private JPanel panel6;
    private JLabel label1;
    private JTextField searchBox;
    private JSplitPane splitPane2;
    private JSplitPane splitPane1;
    private JPanel panel3;
    private JLabel label3;
    private JScrollPane scrollPane1;
    private JList groups;
    private JPanel panel7;
    private JideButton addGroupButton;
    private JPanel panel4;
    private JLabel listLabel;
    private JScrollPane scrollPane2;
    private JList lists;
    private JPanel panel8;
    private JideButton jideButton2;
    private JideButton editButton;
    private JPanel panel5;
    private JLabel label4;
    private JScrollPane scrollPane3;
    private JList loci;
    private JPanel panel9;
    private JideButton jideButton4;
    private JPanel buttonBar;
    private JButton okButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables


    public static void main(String[] args) {
        (new GeneListManagerUI(null)).setVisible(true);
    }
}
