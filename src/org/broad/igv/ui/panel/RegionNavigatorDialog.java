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
 * Created by JFormDesigner on Tue Nov 16 14:34:59 PST 2010
 */

package org.broad.igv.ui.panel;

import org.apache.log4j.Logger;
import org.broad.igv.feature.RegionOfInterest;

import java.awt.*;
import java.awt.event.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Collection;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;

/**
 * @author Damon May
 *         <p/>
 *         This dialog displays a list of RegionOfInterest in a table and allows editing.
 *         Navigation to the start of a region is done by selecting the appropriate row.
 *         <p/>
 *         This dialog is not intended to be persistent.  To view one of these, create it.
 */
public class RegionNavigatorDialog extends JDialog {

    private static Logger log = Logger.getLogger(AttributePanel.class);

    //Column indexes, in case table structure changes
    private static final int TABLE_COLINDEX_CHR = 0;
    private static final int TABLE_COLINDEX_START = 1;
    private static final int TABLE_COLINDEX_END = 2;
    private static final int TABLE_COLINDEX_DESC = 3;


    private DefaultTableModel regionTableModel;
    private List<RegionOfInterest> regions;

    private TableRowSorter<TableModel> regionTableRowSorter;

    public RegionNavigatorDialog(Frame owner) {
        super(owner);
        initComponents();
    }

    public RegionNavigatorDialog(Dialog owner) {
        super(owner);
        initComponents();
    }

    public RegionNavigatorDialog(Frame owner, Collection<RegionOfInterest> regions) {
        this(owner);
        postInit(regions);
    }

    /**
     * Populate the table with the loaded regions
     *
     * @param regionsCollection
     */
    private void postInit(Collection<RegionOfInterest> regionsCollection) {
        this.regions = new ArrayList<RegionOfInterest>(regionsCollection);
        regionTableModel = (DefaultTableModel) regionTable.getModel();

        regionTableModel.setRowCount(regions.size());
        for (int i = 0; i < regions.size(); i++) {
            RegionOfInterest region = regions.get(i);

            regionTableModel.setValueAt(region.getDescription(), i, TABLE_COLINDEX_DESC);
            regionTableModel.setValueAt(region.getDisplayStart(), i, TABLE_COLINDEX_START);
            regionTableModel.setValueAt(region.getDisplayEnd(), i, TABLE_COLINDEX_END);
            regionTableModel.setValueAt(region.getChr(), i, TABLE_COLINDEX_CHR);
        }


        //resize window if small number of regions.  By default, tables are initialized with 20
        //rows.  Only bother resizing if we've got less than 15.  This is hacky and makes the table
        //slightly too tall
        if (regions.size() < 15) {
            int newTableHeight = regions.size() * regionTable.getRowHeight();
            //"extraHeight" is all the height contributed by the top and bottom regions.  There's probably
            //a better way to calculate this.
            //Table initially is 21 rows tall, including 1 row for header
            int extraHeight = getHeight() - (21 * regionTable.getRowHeight());

            int newDialogHeight = newTableHeight + extraHeight;
            if (newDialogHeight < getHeight()) {
                regionTable.setPreferredScrollableViewportSize(new Dimension(regionTable.getPreferredSize().width,
                        newTableHeight));
                setSize(getWidth(), newTableHeight + extraHeight);
                update(getGraphics());
            }
        }


        regionTable.getSelectionModel().addListSelectionListener(new RegionTableSelectionListener());
        regionTableModel.addTableModelListener(new RegionTableModelListener());

        //custom row sorter required for displaying only a subset of rows
        regionTableRowSorter = new TableRowSorter<TableModel>(regionTableModel);
        regionTable.setRowSorter(regionTableRowSorter);

        textFieldSearch.getDocument().addDocumentListener(new SearchFieldDocumentListener());
    }

    private class SearchFieldDocumentListener implements DocumentListener {
        public void changedUpdate(DocumentEvent e) {
            System.err.println("This should not happen");
        }

        public void insertUpdate(DocumentEvent e) {
            showSearchedRegions();
        }

        public void removeUpdate(DocumentEvent e) {
            showSearchedRegions();
        }
    }

    /**
     * Display only those regions whose descriptions match the search string entered by the user
     */
    private void showSearchedRegions() {
        regionTableRowSorter.setRowFilter(new RegionRowFilter(textFieldSearch.getText()));
    }

    /**
     * A row filter that shows only rows that contain filterString, case-insensitive
     */
    private class RegionRowFilter extends RowFilter<TableModel, Object> {
        protected String filterStringLowercase;

        public RegionRowFilter(String filterString) {
            super();
            this.filterStringLowercase = filterString.toLowerCase();
        }

        public boolean include(RowFilter.Entry entry) {
            String regionDesc = (String) entry.getValue(TABLE_COLINDEX_DESC);

            if (filterStringLowercase == null || filterStringLowercase.isEmpty() ||
                    (regionDesc != null &&
                            entry.getValue(TABLE_COLINDEX_DESC).toString().contains(filterStringLowercase)))
                return true;
            return false;
        }
    }

    /**
     * Listen for updates to the cells, save changes to the Regions
     */
    private class RegionTableModelListener implements TableModelListener {
        public void tableChanged(TableModelEvent e) {
            int rowIdx = e.getFirstRow();
            RegionOfInterest region = regions.get(rowIdx);
            switch (e.getColumn()) {
                case TABLE_COLINDEX_DESC:
                    region.setDescription(regionTableModel.getValueAt(rowIdx, TABLE_COLINDEX_DESC).toString());
                    break;
                case TABLE_COLINDEX_START:
                    //stored values are 0-based, viewed values are 1-based.  Check for negative number just in case
                    int storeStartValue =
                            Math.max(0, (Integer) regionTableModel.getValueAt(rowIdx, TABLE_COLINDEX_START) - 1);
                    region.setStart(storeStartValue);
                    break;
                case TABLE_COLINDEX_END:
                    //stored values are 0-based, viewed values are 1-based.  Check for negative number just in case
                    int storeEndValue =
                            Math.max(0, (Integer) regionTableModel.getValueAt(rowIdx, TABLE_COLINDEX_END) - 1);
                    region.setEnd(storeEndValue);
                    break;
            }

        }
    }

    /**
     * Listen for selection change events, navigate UI to start of selected region
     */
    private class RegionTableSelectionListener implements ListSelectionListener {
        public void valueChanged(ListSelectionEvent e) {
            if (!e.getValueIsAdjusting()) {
                int[] selectedRows = regionTable.getSelectedRows();
                if (selectedRows != null && selectedRows.length == 1) {
                    int selectedModelRow = regionTableRowSorter.convertRowIndexToModel(selectedRows[0]);
                    RegionOfInterest region = regions.get(selectedModelRow);
                    log.debug("Nav to region " + region.getDescription() + "chr=" + region.getChr() +
                            ", start=" + region.getStart() + ", end=" + region.getEnd());

                    // Option (1), center on the region without changing resolution
                    //FrameManager.getDefaultFrame().centerOnLocation(region.getChr(), region.getCenter());

                    // Option (2), zoom and center on region, with an interval equal to the region length on
                    // either side for context
                    int start = region.getStart() - region.getLength();
                    int end = region.getEnd() + region.getLength();
                    FrameManager.getDefaultFrame().jumpTo(region.getChr(), start, end);

                }
            }
        }
    }


    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        scrollPane1 = new JScrollPane();
        regionTable = new JTable();
        panel1 = new JPanel();
        textFieldSearch = new JTextField();
        label1 = new JLabel();
        button1 = new JButton();
        cancelAction = new CancelAction();

        //======== this ========
        setTitle("Regions of Interest");
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(null);
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(new BorderLayout());

                //======== scrollPane1 ========
                {

                    //---- regionTable ----
                    regionTable.setModel(new DefaultTableModel(
                            new Object[][]{
                                    {null, null, null, null},
                            },
                            new String[]{
                                    "Chr", "Start", "End", "Description"
                            }
                    ) {
                        Class<?>[] columnTypes = new Class<?>[]{
                                String.class, Integer.class, Integer.class, Object.class
                        };
                        boolean[] columnEditable = new boolean[]{
                                false, true, true, true
                        };

                        @Override
                        public Class<?> getColumnClass(int columnIndex) {
                            return columnTypes[columnIndex];
                        }

                        @Override
                        public boolean isCellEditable(int rowIndex, int columnIndex) {
                            return columnEditable[columnIndex];
                        }
                    });
                    {
                        TableColumnModel cm = regionTable.getColumnModel();
                        cm.getColumn(0).setPreferredWidth(50);
                        cm.getColumn(1).setPreferredWidth(100);
                        cm.getColumn(2).setPreferredWidth(100);
                        cm.getColumn(3).setPreferredWidth(200);
                    }
                    regionTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
                    regionTable.setAutoCreateRowSorter(true);
                    scrollPane1.setViewportView(regionTable);
                }
                contentPanel.add(scrollPane1, BorderLayout.CENTER);

                //======== panel1 ========
                {
                    panel1.setLayout(new BorderLayout());
                    panel1.add(textFieldSearch, BorderLayout.CENTER);

                    //---- label1 ----
                    label1.setText("Search");
                    panel1.add(label1, BorderLayout.WEST);

                    //---- button1 ----
                    button1.setAction(cancelAction);
                    button1.setText("Clear");
                    panel1.add(button1, BorderLayout.EAST);
                }
                contentPanel.add(panel1, BorderLayout.SOUTH);
            }
            dialogPane.add(contentPanel, BorderLayout.CENTER);
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
    private JTable regionTable;
    private JPanel panel1;
    private JTextField textFieldSearch;
    private JLabel label1;
    private JButton button1;
    private CancelAction cancelAction;
    // JFormDesigner - End of variables declaration  //GEN-END:variables


    private class SaveAction extends AbstractAction {
        private SaveAction() {
            // JFormDesigner - Action initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
            // Generated using JFormDesigner non-commercial license
            // JFormDesigner - End of action initialization  //GEN-END:initComponents
        }

        public void actionPerformed(ActionEvent e) {
            // TODO add your code here
        }
    }

    private class CancelAction extends AbstractAction {
        private CancelAction() {
            // JFormDesigner - Action initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
            // Generated using JFormDesigner non-commercial license
            putValue(NAME, "Cancel");
            putValue(SHORT_DESCRIPTION, "Clear search box");
            // JFormDesigner - End of action initialization  //GEN-END:initComponents
        }

        public void actionPerformed(ActionEvent e) {
            textFieldSearch.setText("");
        }
    }
}
