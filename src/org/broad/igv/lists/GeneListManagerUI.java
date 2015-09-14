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
 * Created by JFormDesigner on Mon Dec 13 14:39:17 EST 2010
 */

package org.broad.igv.lists;

import javax.swing.plaf.*;
import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.cbio.FilterGeneNetworkUI;
import org.broad.igv.gitools.Gitools;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.igv.ui.util.MessageUtils;

import javax.swing.*;
import javax.swing.border.BevelBorder;
import javax.swing.border.EmptyBorder;
import javax.swing.border.EtchedBorder;
import javax.swing.border.LineBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.List;

/**
 * @author Jim RObinson
 *         <p/>
 *         public String getToolTipText(MouseEvent evt) {
 *         // Get item index
 *         int index = locationToIndex(evt.getPoint());
 *         // Get item Object i
 *         tem = getModel().getElementAt(index);
 *         // Return the tool tip text
 *         return "tool tip for "+item; }
 */
public class GeneListManagerUI extends JDialog {

    private static Logger log = Logger.getLogger(GeneListManagerUI.class);

    private static String ALL = "All";
    private static String DEFAULT_ACTION_TEXT = "View";
    private static final GeneListListener DEFAULT_ACTION_LISTENER = new GeneListListener() {
        @Override
        public void actionPerformed(JDialog dialog, GeneList geneList) {
            IGV.getInstance().setGeneList(geneList);
            dialog.setVisible(false);
            dialog.dispose();
        }
    };

    // This is a static so its "remembered" withing the context of a session.
    static String lastSelectedGroup;

    String selectedGroup;
    ListListModel listModel;
    GeneListModel geneListModel;

    /**
     * Map of gene list name -> gene list
     */
    Map<String, GeneList> geneLists;
    private String selectedList;

    GeneListManager manager;

    private static GeneListManagerUI instance;
    private GeneListListener listener;

    /**
     * Get or create new GeneListManagerUI.
     * If old instance retrieved, the owner IS NOT changed
     *
     * @param owner
     * @return
     */
    public static GeneListManagerUI getInstance(Frame owner) {
        return getInstance(owner, DEFAULT_ACTION_TEXT, DEFAULT_ACTION_LISTENER);
    }

    /**
     * @param owner
     * @param actionText Text to be displayed when the action is performed. Default is "View", for simply viewing the gene list
     * @param listener   The action performed with the specified GeneList
     * @return
     */
    public static GeneListManagerUI getInstance(Frame owner, String actionText, GeneListListener listener) {
        if (instance == null) {
            instance = new GeneListManagerUI(owner);
        }

        instance.actionButton.setText(actionText);
        instance.listener = listener;

        if (listener != DEFAULT_ACTION_LISTENER) {
            instance.viewNetworkButton.setVisible(false);
        }

        return instance;
    }


    private GeneListManagerUI(Frame owner) {
        super(owner);
        manager = GeneListManager.getInstance();
        initComponents();
        initLists();

//        boolean showTDMButton = Boolean.parseBoolean(System.getProperty("enableGitools", "false"));
//        exportTDMButton.setVisible(showTDMButton);
    }

    private void initLists() {
        geneLists = manager.getGeneLists();

        groupJList.setModel(new AbstractListModel() {

            ArrayList<String> groups = new ArrayList();

            {
                groups.add(ALL);
                groups.addAll(manager.getGroups());
            }


            public int getSize() {
                return groups.size();
            }

            public Object getElementAt(int i) {
                return groups.get(i);
            }
        });

        listModel = new ListListModel();
        glJList.setModel(listModel);

        geneListModel = new GeneListModel();
        lociJList.setModel(geneListModel);

        if (lastSelectedGroup != null) {
            groupJList.setSelectedValue(lastSelectedGroup, true);
        } else {
            groupJList.setSelectedIndex(0);
        }
        glJList.setSelectedIndex(0);
    }


    private void groupsValueChanged(ListSelectionEvent e) {
        selectedGroup = (String) groupJList.getSelectedValue();
        lastSelectedGroup = selectedGroup;
        deleteGroupButton.setEnabled(!(GeneListManager.USER_GROUP.equals(selectedGroup) ||
                !GeneListManager.DEFAULT_GENE_LISTS.contains(selectedGroup)));
        updateListModel();
    }

    private void listsValueChanged(ListSelectionEvent e) {

        selectedList = (String) glJList.getSelectedValue();
        if (selectedList == null) {
            geneListModel.clear();
        } else {
            actionButton.setEnabled(true);
            GeneList gl = listModel.getGeneList(selectedList);
            geneListModel.setGeneList(gl);
            lociJList.setModel(geneListModel);
            lociJList.updateUI();

            editButton.setEnabled(gl.isEditable());
            //deleteButton.setEnabled(gl.isEditable());
        }
    }

    private void listLabelMouseClicked(MouseEvent e) {
        listModel.sort();
        glJList.updateUI();
    }

    private void searchBoxKeyReleased(KeyEvent e) {

        updateListModel();
    }

    private void updateListModel() {
        listModel.filter();
        glJList.clearSelection();
        glJList.updateUI();
        lociJList.updateUI();
    }

    /**
     * Stub method generated by JFormDesigner t0 place custom code
     */
    private void createUIComponents() {
        glJList = new JList() {
            public String getToolTipText(MouseEvent evt) {
                int index = locationToIndex(evt.getPoint());
                if (index >= 0) {
                    Object item = getModel().getElementAt(index);
                    GeneList gl = geneLists.get(item);
                    return gl.getDescription();
                }
                return null;
            }
        };
    }


    private void editButtonActionPerformed(ActionEvent e) {

        String selection = (String) glJList.getSelectedValue();
        if (selection != null) {
            GeneList geneList = geneLists.get(selection);
            GeneListEditDialog dlg = new GeneListEditDialog(this, geneList);
            dlg.setVisible(true);

            if (!dlg.isCanceled()) {
                geneListModel.setGeneList(geneList);
                listModel.filter();
                glJList.updateUI();
                lociJList.updateUI();
            }

        }
    }


    private void newListActionPerformed(ActionEvent e) {
        GeneList geneList = new GeneList();
        GeneListEditDialog dlg = new GeneListEditDialog(this, geneList);
        dlg.setVisible(true);

        if (!dlg.isCanceled()) {
            geneList.setGroup(manager.USER_GROUP);
            manager.addGeneList(geneList);
            listModel.add(geneList);
            glJList.updateUI();
            glJList.setSelectedValue(geneList.getName(), true);
            groupJList.updateUI();
            //loci.updateUI();
        }
    }


    private void copyListButtonActionPerformed(ActionEvent e) {
        String selection = (String) glJList.getSelectedValue();
        if (selection != null) {
            GeneList geneList = geneLists.get(selection);
            GeneList copiedList = geneList.copy();
            GeneListEditDialog dlg = new GeneListEditDialog(this, copiedList);
            dlg.setVisible(true);

            if (!dlg.isCanceled()) {
                listModel.add(copiedList);
                glJList.updateUI();
                glJList.setSelectedValue(copiedList.getName(), true);
                //loci.updateUI();
            }
        }
    }

    private void deleteButtonActionPerformed(ActionEvent e) {
        String selection = (String) glJList.getSelectedValue();
        if (selection != null) {
            if (MessageUtils.confirm("<html>Are you sure you want to delete list '" + selection + "' ?" +
                    "<br>associated lists? &nbsp;<b>This action cannot be undone")) {

                boolean groupRemoved = manager.deleteList(selection);
                if (groupRemoved) {
                    initLists();
                } else {
                    geneLists.remove(selection);
                    listModel.filter();
                    glJList.updateUI();
                    glJList.setSelectedIndex(0);
                    lociJList.updateUI();
                }
            }
        }
    }

    /**
     * Import a gene list group from a GMT file
     *
     * @param e
     */
    private void importButtonActionPerformed(ActionEvent e) {

        File gmtFile = FileDialogUtils.chooseFile("Import GMT File");
        if (gmtFile != null) {

            try {
                List<GeneList> loadedLists = manager.importFile(gmtFile);
                initLists();

                if (loadedLists.size() > 0) {
                    groupJList.setSelectedValue(loadedLists.get(0).getGroup(), true);
                    glJList.setSelectedValue(loadedLists.get(0).getName(), true);
                    groupJList.updateUI();
                    glJList.updateUI();
                }

            } catch (IOException e1) {
                log.error("Error updating genome property file", e1);
                MessageUtils.showMessage("Error importing .gmt file: " + e1.getMessage());
            }
        }
    }

    /**
     * Export the selected group to a GMT file
     *
     * @param e
     */
    private void exportButtonActionPerformed(ActionEvent e) {
        if (selectedGroup != null) {
            File userDir = DirectoryManager.getUserDirectory();
            File initFile = new File(selectedGroup + ".gmt");
            File glFile = FileDialogUtils.chooseFile("Save gene lists", userDir, initFile, FileDialogUtils.SAVE);
            if (glFile != null) {
                try {
                    manager.exportGMT(selectedGroup, glFile);
                } catch (Exception e1) {
                    log.error("Error exporting gene lists", e1);
                    MessageUtils.showMessage("Error exporting gene lists: " + e1.getMessage());
                }
            }
        }
    }


    private void deleteGroupButtonActionPerformed(ActionEvent e) {

        if (selectedGroup != null && MessageUtils.confirm(this, "<html>Are you sure you want to delete group '" +
                selectedGroup + "' and all " +
                "<br>associated lists? &nbsp;<b>This action cannot be undone")) {
            manager.deleteGroup(selectedGroup);
            initLists();
        }
    }


    private void closeButtonActionPerformed(ActionEvent e) {
        setVisible(false);
    }

    private void actionButtonActionPerformed(ActionEvent e) {
        if (selectedList != null) {
            GeneList geneList = geneLists.get(selectedList);
            listener.actionPerformed(this, geneList);
        }
    }

    private void glJListMouseClicked(MouseEvent e) {
        if (e.getClickCount() > 1) {
            this.actionButtonActionPerformed(null);
        }
    }

    private void retrieveNetworkButtonActionPerformed(ActionEvent e) {
        if (selectedList != null) {
            GeneList geneList = geneLists.get(selectedList);
            FilterGeneNetworkUI fgnUI = new FilterGeneNetworkUI(IGV.getMainFrame(), geneList);
            fgnUI.setVisible(true);
        }
    }

    private void exportTDMButtonActionPerformed_deleteme(ActionEvent e) {
        if (selectedList != null) {
            GeneList geneList = geneLists.get(selectedList);

            File file = FileDialogUtils.chooseFile("Export TDM file", null, FileDialogUtils.SAVE);
            if (file != null) {
                try {
                    Gitools.exportTDM(geneList.getLoci(), file);
                } catch (IOException exc) {
                    log.error("Error exporting TDM data", exc);
                    MessageUtils.showErrorMessage("Error exporting TDM data", exc);
                }
            }
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

            if (selectedGroup != null && !selectedGroup.equals(ALL)) {
                if (!geneList.getGroup().equals(selectedGroup)) {
                    return false;
                }
            }

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

        public void add(GeneList geneList) {
            geneLists.put(geneList.getName(), geneList);
            if (isPassFilter(geneList)) {
                filteredNames.add(geneList.getName());
            }
        }
    }

    class GeneListModel extends AbstractListModel {

        java.util.List<String> genes;

        GeneListModel() {
            genes = new ArrayList();
        }

        void setGeneList(GeneList geneList) {
            genes = geneList == null ? new ArrayList() : new ArrayList(geneList.getLoci());
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
        createUIComponents();

        dialogPane = new JPanel();
        contentPanel = new JPanel();
        panel2 = new JPanel();
        panel1 = new JPanel();
        panel6 = new JPanel();
        label1 = new JLabel();
        searchBox = new JTextField();
        panel10 = new JPanel();
        splitPane2 = new JSplitPane();
        splitPane1 = new JSplitPane();
        panel3 = new JPanel();
        label3 = new JLabel();
        scrollPane1 = new JScrollPane();
        groupJList = new JList();
        panel7 = new JPanel();
        importButton = new JButton();
        exportButton = new JButton();
        deleteGroupButton = new JButton();
        panel4 = new JPanel();
        listLabel = new JLabel();
        scrollPane2 = new JScrollPane();
        panel8 = new JPanel();
        newList = new JButton();
        copyListButton = new JButton();
        editButton = new JButton();
        deleteButton = new JButton();
        panel5 = new JPanel();
        label4 = new JLabel();
        scrollPane3 = new JScrollPane();
        lociJList = new JList();
        panel9 = new JPanel();
        buttonBar = new JPanel();
        exportTDMButton = new JButton();
        viewNetworkButton = new JButton();
        actionButton = new JButton();
        closeButton = new JButton();

        //======== this ========
        setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
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
                            panel6.setMinimumSize(new Dimension(0, 0));
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

                        //======== panel10 ========
                        {
                            panel10.setPreferredSize(new Dimension(5, 0));
                            panel10.setMinimumSize(new Dimension(5, 0));
                            panel10.setOpaque(false);
                            panel10.setLayout(null);
                        }
                        panel1.add(panel10);
                    }
                    panel2.add(panel1);
                }
                contentPanel.add(panel2, BorderLayout.NORTH);

                //======== splitPane2 ========
                {
                    splitPane2.setBorder(null);
                    splitPane2.setDividerLocation(600);

                    //======== splitPane1 ========
                    {
                        splitPane1.setBorder(null);
                        splitPane1.setDividerLocation(270);

                        //======== panel3 ========
                        {
                            panel3.setBorder(LineBorder.createBlackLineBorder());
                            panel3.setLayout(new BorderLayout());

                            //---- label3 ----
                            label3.setText("Group");
                            panel3.add(label3, BorderLayout.NORTH);

                            //======== scrollPane1 ========
                            {
                                scrollPane1.setPreferredSize(new Dimension(80, 140));

                                //---- groupJList ----
                                groupJList.setModel(new AbstractListModel() {
                                    String[] values = {
                                        "All"
                                    };
                                    @Override
                                    public int getSize() { return values.length; }
                                    @Override
                                    public Object getElementAt(int i) { return values[i]; }
                                });
                                groupJList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
                                groupJList.addListSelectionListener(new ListSelectionListener() {
                                    @Override
                                    public void valueChanged(ListSelectionEvent e) {
                                        groupsValueChanged(e);
                                    }
                                });
                                scrollPane1.setViewportView(groupJList);
                            }
                            panel3.add(scrollPane1, BorderLayout.CENTER);

                            //======== panel7 ========
                            {
                                panel7.setLayout(new BoxLayout(panel7, BoxLayout.X_AXIS));

                                //---- importButton ----
                                importButton.setIcon(null);
                                importButton.setText("Import");
                                importButton.setToolTipText("Import a .gmt file, .grp, or .bed file");
                                importButton.addActionListener(new ActionListener() {
                                    @Override
                                    public void actionPerformed(ActionEvent e) {
                                        importButtonActionPerformed(e);
                                    }
                                });
                                panel7.add(importButton);

                                //---- exportButton ----
                                exportButton.setText("Export");
                                exportButton.addActionListener(new ActionListener() {
                                    @Override
                                    public void actionPerformed(ActionEvent e) {
                                        exportButtonActionPerformed(e);
                                    }
                                });
                                panel7.add(exportButton);

                                //---- deleteGroupButton ----
                                deleteGroupButton.setText("Delete");
                                deleteGroupButton.setEnabled(false);
                                deleteGroupButton.addActionListener(new ActionListener() {
                                    @Override
                                    public void actionPerformed(ActionEvent e) {
                                        deleteGroupButtonActionPerformed(e);
                                    }
                                });
                                panel7.add(deleteGroupButton);
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

                                //---- glJList ----
                                glJList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
                                glJList.addListSelectionListener(new ListSelectionListener() {
                                    @Override
                                    public void valueChanged(ListSelectionEvent e) {
                                        listsValueChanged(e);
                                    }
                                });
                                glJList.addMouseListener(new MouseAdapter() {
                                    @Override
                                    public void mouseClicked(MouseEvent e) {
                                        glJListMouseClicked(e);
                                    }
                                });
                                scrollPane2.setViewportView(glJList);
                            }
                            panel4.add(scrollPane2, BorderLayout.CENTER);

                            //======== panel8 ========
                            {
                                panel8.setLayout(new BoxLayout(panel8, BoxLayout.X_AXIS));

                                //---- newList ----
                                newList.setIcon(null);
                                newList.setText("New");
                                newList.addActionListener(new ActionListener() {
                                    @Override
                                    public void actionPerformed(ActionEvent e) {
                                        newListActionPerformed(e);
                                    }
                                });
                                panel8.add(newList);

                                //---- copyListButton ----
                                copyListButton.setText("Copy");
                                copyListButton.addActionListener(new ActionListener() {
                                    @Override
                                    public void actionPerformed(ActionEvent e) {
                                        copyListButtonActionPerformed(e);
                                    }
                                });
                                panel8.add(copyListButton);

                                //---- editButton ----
                                editButton.setText("Edit");
                                editButton.setEnabled(false);
                                editButton.addActionListener(new ActionListener() {
                                    @Override
                                    public void actionPerformed(ActionEvent e) {
                                        editButtonActionPerformed(e);
                                    }
                                });
                                panel8.add(editButton);

                                //---- deleteButton ----
                                deleteButton.setIcon(null);
                                deleteButton.setText("Delete");
                                deleteButton.addActionListener(new ActionListener() {
                                    @Override
                                    public void actionPerformed(ActionEvent e) {
                                        deleteButtonActionPerformed(e);
                                    }
                                });
                                panel8.add(deleteButton);
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

                            //---- lociJList ----
                            lociJList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
                            lociJList.setSelectionBackground(Color.white);
                            lociJList.setSelectionForeground(Color.black);
                            lociJList.setFocusable(false);
                            scrollPane3.setViewportView(lociJList);
                        }
                        panel5.add(scrollPane3, BorderLayout.CENTER);

                        //======== panel9 ========
                        {
                            panel9.setLayout(new BoxLayout(panel9, BoxLayout.X_AXIS));
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
                buttonBar.setBorder(null);
                buttonBar.setLayout(new FlowLayout(FlowLayout.RIGHT));

                //---- exportTDMButton ----
                exportTDMButton.setText("Export TDM");
                exportTDMButton.setVisible(false);
                exportTDMButton.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        exportTDMButtonActionPerformed_deleteme(e);
                    }
                });
                buttonBar.add(exportTDMButton);

                //---- viewNetworkButton ----
                viewNetworkButton.setText("Retrieve Network");
                viewNetworkButton.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        retrieveNetworkButtonActionPerformed(e);
                    }
                });
                buttonBar.add(viewNetworkButton);

                //---- actionButton ----
                actionButton.setText("View");
                actionButton.setEnabled(false);
                actionButton.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        actionButtonActionPerformed(e);
                    }
                });
                buttonBar.add(actionButton);

                //---- closeButton ----
                closeButton.setText("Close");
                closeButton.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        closeButtonActionPerformed(e);
                    }
                });
                buttonBar.add(closeButton);
            }
            dialogPane.add(buttonBar, BorderLayout.SOUTH);
        }
        contentPane.add(dialogPane, BorderLayout.CENTER);
        setSize(845, 580);
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
    private JPanel panel10;
    private JSplitPane splitPane2;
    private JSplitPane splitPane1;
    private JPanel panel3;
    private JLabel label3;
    private JScrollPane scrollPane1;
    private JList groupJList;
    private JPanel panel7;
    private JButton importButton;
    private JButton exportButton;
    private JButton deleteGroupButton;
    private JPanel panel4;
    private JLabel listLabel;
    private JScrollPane scrollPane2;
    private JList glJList;
    private JPanel panel8;
    private JButton newList;
    private JButton copyListButton;
    private JButton editButton;
    private JButton deleteButton;
    private JPanel panel5;
    private JLabel label4;
    private JScrollPane scrollPane3;
    private JList lociJList;
    private JPanel panel9;
    private JPanel buttonBar;
    private JButton exportTDMButton;
    private JButton viewNetworkButton;
    private JButton actionButton;
    private JButton closeButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables

    /**
     * Custom action listener for receiving a GeneList
     * and taking some action based on that. Listeners
     * implementing this interface can be used as callbacks
     * once a GeneList is selected from GeneListManagerUI
     *
     * @author jacob
     * @date 2013-Jan-04
     * @api
     */
    public interface GeneListListener {

        void actionPerformed(JDialog dialog, GeneList geneList);
    }
}