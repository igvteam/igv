/**
 * Copyright (c) 2010-2011 by Fred Hutchinson Cancer Research Center.  All Rights Reserved.

 * This software is licensed under the terms of the GNU Lesser General
 * Public License (LGPL), Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.

 * THE SOFTWARE IS PROVIDED "AS IS." FRED HUTCHINSON CANCER RESEARCH CENTER MAKES NO
 * REPRESENTATIONS OR WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED,
 * INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 * PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS,
 * WHETHER OR NOT DISCOVERABLE.  IN NO EVENT SHALL FRED HUTCHINSON CANCER RESEARCH
 * CENTER OR ITS TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR
 * ANY DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR
 * CONSEQUENTIAL DAMAGES, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS,
 * REGARDLESS OF  WHETHER FRED HUTCHINSON CANCER RESEARCH CENTER SHALL BE ADVISED,
 * SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE
 * FOREGOING.
 */
package org.broad.igv.ui.panel;

import org.apache.log4j.Logger;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.feature.SequenceManager;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.IGVMainFrame;
import org.broad.igv.ui.util.MessageUtils;

import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.*;
import java.util.*;
import java.util.List;
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
 *
 * todo: tell the regions observable that this object is toast, when it goes away
 */
public class RegionNavigatorDialog extends JDialog implements Observer{

    private static Logger log = Logger.getLogger(AttributePanel.class);

    //Column indexes, in case table structure changes
    private static final int TABLE_COLINDEX_CHR = 0;
    private static final int TABLE_COLINDEX_START = 1;
    private static final int TABLE_COLINDEX_END = 2;
    private static final int TABLE_COLINDEX_DESC = 3;

    //The active instance of RegionNavigatorDialog (only one at a time)
    public static RegionNavigatorDialog activeInstance;


    private DefaultTableModel regionTableModel;
//    private List<RegionOfInterest> regions;

    private TableRowSorter<TableModel> regionTableRowSorter;

    //Indicates that we're in the process of synching the table with the regions list, so we shouldn't
    //do anything about TableChanged events.
    private boolean synchingRegions = false;

    /**
     * Return the active RegionNavigatorDialog, or null if none.
     * @return
     */
    public static RegionNavigatorDialog getActiveInstance() {
        return activeInstance;
    }

    /**
     * dispose the active instance and get rid of the pointer. Return whether or not there was an
     * active instance
     */
    public static boolean destroyActiveInstance() {
        if (activeInstance == null)
            return false;
        activeInstance.dispose();
        activeInstance = null;
        return true;
    }

    public RegionNavigatorDialog(Frame owner) {
        super(owner);
        initComponents();
        postInit();
    }

    public RegionNavigatorDialog(Dialog owner) {
        super(owner);
        initComponents();
        postInit();
    }

    public void update(Observable observable, Object object)
    {
        synchRegions();
    }

    /**
     * Synchronize the regions ArrayList with the passed-in regionsCollection, and update UI
     */
    public void synchRegions()
    {
        //Indicate that we're synching regions, so that we don't respond to tableChanged events
        synchingRegions = true;
        List<RegionOfInterest> regions = retrieveRegionsAsList();
        regionTableModel = (DefaultTableModel) regionTable.getModel();
        while (regionTableModel.getRowCount() > 0)
            regionTableModel.removeRow(0);
        regionTableModel.setRowCount(regions.size());
        for (int i = 0; i < regions.size(); i++) {
            RegionOfInterest region = regions.get(i);

            regionTableModel.setValueAt(region.getDescription(), i, TABLE_COLINDEX_DESC);
            regionTableModel.setValueAt(region.getDisplayStart(), i, TABLE_COLINDEX_START);
            regionTableModel.setValueAt(region.getDisplayEnd(), i, TABLE_COLINDEX_END);
            regionTableModel.setValueAt(region.getChr(), i, TABLE_COLINDEX_CHR);
        }
        //Done synching regions, allow ourselves to respond to tableChanged events
        synchingRegions = false;

        regionTableModel.fireTableDataChanged();
    }

    private List<RegionOfInterest> retrieveRegionsAsList()
    {
        return new ArrayList<RegionOfInterest>(IGVMainFrame.getInstance().getSession().getAllRegionsOfInterest());
    }

    /**
     * Populate the table with the loaded regions
     *
     */
    private void postInit() {
        regionTableModel = (DefaultTableModel) regionTable.getModel();

        regionTable.getSelectionModel().addListSelectionListener(new RegionTableSelectionListener());
        regionTableModel.addTableModelListener(new RegionTableModelListener());

        //custom row sorter required for displaying only a subset of rows
        regionTableRowSorter = new TableRowSorter<TableModel>(regionTableModel);
        regionTable.setRowSorter(regionTableRowSorter);
        regionTableRowSorter.setRowFilter(new RegionRowFilter());

        textFieldSearch.getDocument().addDocumentListener(new SearchFieldDocumentListener());

        activeInstance = this;
        updateChromosomeDisplayed();

        synchRegions();

        IGVMainFrame.getInstance().getSession().getRegionsOfInterestObservable().addObserver(this);

        //resize window if small number of regions.  By default, tables are initialized with 20
        //rows, and that can look ungainly for empty windows or windows with a few rows.
        //This correction is rather hacky. Minimum size of 5 rows set.
        int newTableHeight = Math.min(regionTableModel.getRowCount()+1,5) * regionTable.getRowHeight();
        //This is quite hacky -- need to find the size of the other components programmatically somehow, since
        //it will vary on different platforms
        int extraHeight = 225;

        int newDialogHeight = newTableHeight + extraHeight;
        if (newDialogHeight < getHeight()) {
            regionTable.setPreferredScrollableViewportSize(new Dimension(regionTable.getPreferredSize().width,
                    newTableHeight));
            setSize(getWidth(), newTableHeight + extraHeight);
            update(getGraphics());
        }

        regionTable.addMouseListener(new RegionTablePopupHandler());
    }

    private class SearchFieldDocumentListener implements DocumentListener {
        public void changedUpdate(DocumentEvent e) {
            System.err.println("This should not happen");
        }

        public void insertUpdate(DocumentEvent e) {
            regionTableModel.fireTableDataChanged();
        }

        public void removeUpdate(DocumentEvent e) {
            regionTableModel.fireTableDataChanged();
        }
    }

    /**
     * When chromosome that's displayed is changed, need to update displayed regions.  showSearchedRegions will do that
     */
    public void updateChromosomeDisplayed()
    {
//        regionTable.updateUI();
//        showSearchedRegions();
        regionTableModel.fireTableDataChanged();
    }

    /**
     * Test whether we should display an entry
     * @param regionChr
     * @param regionDesc
     * @return
     */
    protected boolean shouldIncludeRegion(String regionChr, String regionDesc)
    {
            //if table is empty, a non-region event is fed here.  Test for it and don't display
            if (regionChr == null)
                return false;

            String filterStringLowercase = null;
            if (textFieldSearch.getText() != null)
                filterStringLowercase = textFieldSearch.getText().toLowerCase();

            String chr = FrameManager.getDefaultFrame().getChrName();

            //show only regions matching the search string (if specified)
            if ((filterStringLowercase != null && !filterStringLowercase.isEmpty() &&
                (regionDesc == null || !regionDesc.toLowerCase().contains(filterStringLowercase))))
                return false;

            //if this checkbox is checked, show all chromosomes
            if (checkBoxShowAllChrs.isSelected())
                return true;

            //show only regions in the current chromosome
            if (chr != null && !chr.isEmpty() && !regionChr.equals(chr))
                return false;


            return true;
    }

    /**
     * A row filter that shows only rows that contain filterString, case-insensitive
     */
    private class RegionRowFilter extends RowFilter<TableModel, Object> {

        public RegionRowFilter() {
            super();
        }

        public boolean include(RowFilter.Entry entry) {
            return shouldIncludeRegion((String) entry.getValue(TABLE_COLINDEX_CHR),
                    (String) entry.getValue(TABLE_COLINDEX_DESC));
        }
    }

    /**
     * Listen for updates to the cells, save changes to the Regions
     */
    private class RegionTableModelListener implements TableModelListener {
        public void tableChanged(TableModelEvent e) {
            //If we're in the middle of synching regions, do nothing
            if (synchingRegions)
                return;

            List<RegionOfInterest> regions = retrieveRegionsAsList();
            int firstRow = e.getFirstRow();
            //range checking because this method gets called after a clear event, and we don't want to
            //try to find an updated region then
            if (e.getFirstRow() > regions.size()-1)
                return;

            //must convert row index from view to model, in case of sorting, filtering
            int rowIdx = 0;

            try
            {
                rowIdx = regionTable.getRowSorter().convertRowIndexToModel(e.getFirstRow());
            }
            catch (ArrayIndexOutOfBoundsException x)
            {
                return;
            }

            RegionOfInterest region = regions.get(rowIdx);
            switch (e.getColumn()) {
                case TABLE_COLINDEX_DESC:
                    Object descObject = regionTableModel.getValueAt(rowIdx, TABLE_COLINDEX_DESC);
                    if (descObject != null)
                         region.setDescription(descObject.toString());
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
                List<RegionOfInterest> regions = retrieveRegionsAsList();

                int[] selectedRows = regionTable.getSelectedRows();
                if (selectedRows != null && selectedRows.length > 0
                        && regions.size() >= selectedRows.length)  //dhmay: this is hacky. Bad things can happen with clear regions
                {
                    RegionOfInterest firstStartRegion = null;
                    RegionOfInterest lastEndRegion = null;

                    Set<String> selectedChrs = new HashSet<String>();

                    //Figure out which region has the first start and which has the last end
                    for (int selectedRowIndex : selectedRows)
                    {
                        int selectedModelRow = regionTableRowSorter.convertRowIndexToModel(selectedRowIndex);
                        RegionOfInterest region = regions.get(selectedModelRow);
                        selectedChrs.add(region.getChr());
                        if (firstStartRegion == null || region.getStart() < firstStartRegion.getStart())
                            firstStartRegion = region;
                        if (lastEndRegion == null || region.getEnd() > lastEndRegion.getEnd())
                            lastEndRegion = region;
                    }

                    //If there are multiple chromosomes represented in the selection, do nothing.
                    //Because what would we do? Maybe a status message should be displayed somehow, but a
                    //dialog would get annoying.
                    if (selectedChrs.size() > 1)
                        return;

                    if (checkBoxZoomWhenNav.isSelected())
                    {
                        // Option (1), zoom and center on group of selected regions, with an interval equal to
                        // 20% of the length of the end regions on either side for context (dhmay reduced from 100%)
                        int start = firstStartRegion.getStart() - (int) (0.2 * firstStartRegion.getLength());
                        int end = lastEndRegion.getEnd() + (int) (0.2 * lastEndRegion.getLength());
                        FrameManager.getDefaultFrame().jumpTo(selectedChrs.iterator().next(), start, end);
                    }
                    else
                    {
                        // Option (2), center on the FIRST selected region without changing resolution
                        FrameManager.getDefaultFrame().centerOnLocation(firstStartRegion.getCenter());
                    }

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
        panel2 = new JPanel();
        buttonAddRegion = new JButton();
        button2 = new JButton();
        panel3 = new JPanel();
        checkBoxZoomWhenNav = new JCheckBox();
        checkBoxShowAllChrs = new JCheckBox();
        cancelAction = new CancelAction();
        addRegionAction = new AddRegionAction();
        actionRemoveRegions = new RemoveSelectedRegionsAction();
        showAllChromosomesCheckboxAction = new ShowAllChromosomesCheckboxAction();

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
                        new Object[][] {
                            {null, null, null, null},
                        },
                        new String[] {
                            "Chr", "Start", "End", "Description"
                        }
                    ) {
                        Class<?>[] columnTypes = new Class<?>[] {
                            String.class, Integer.class, Integer.class, Object.class
                        };
                        boolean[] columnEditable = new boolean[] {
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
                    regionTable.setAutoCreateRowSorter(true);
                    scrollPane1.setViewportView(regionTable);
                }
                contentPanel.add(scrollPane1, BorderLayout.CENTER);

                //======== panel1 ========
                {
                    panel1.setLayout(new BorderLayout());

                    //---- textFieldSearch ----
                    textFieldSearch.setToolTipText("When navigating to a region, zoom to the region?");
                    panel1.add(textFieldSearch, BorderLayout.CENTER);

                    //---- label1 ----
                    label1.setText("Search");
                    panel1.add(label1, BorderLayout.WEST);

                    //---- button1 ----
                    button1.setAction(cancelAction);
                    button1.setText("Clear Search");
                    panel1.add(button1, BorderLayout.EAST);

                    //======== panel2 ========
                    {
                        panel2.setLayout(new BorderLayout());

                        //---- buttonAddRegion ----
                        buttonAddRegion.setAction(addRegionAction);
                        buttonAddRegion.setText("Add");
                        panel2.add(buttonAddRegion, BorderLayout.WEST);

                        //---- button2 ----
                        button2.setAction(actionRemoveRegions);
                        button2.setText("Remove Selected");
                        panel2.add(button2, BorderLayout.CENTER);

                        //======== panel3 ========
                        {
                            panel3.setLayout(new BorderLayout());

                            //---- checkBoxZoomWhenNav ----
                            checkBoxZoomWhenNav.setText("Zoom to Region");
                            checkBoxZoomWhenNav.setToolTipText("When navigating to a region, change zoom level?");
                            checkBoxZoomWhenNav.setSelected(true);
                            panel3.add(checkBoxZoomWhenNav, BorderLayout.CENTER);

                            //---- checkBoxShowAllChrs ----
                            checkBoxShowAllChrs.setAction(showAllChromosomesCheckboxAction);
                            checkBoxShowAllChrs.setToolTipText("View regions from all chromosomes (othrwise, current chromosome only)");
                            checkBoxShowAllChrs.setSelected(true);
                            panel3.add(checkBoxShowAllChrs, BorderLayout.WEST);
                        }
                        panel2.add(panel3, BorderLayout.EAST);
                    }
                    panel1.add(panel2, BorderLayout.NORTH);
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
    private JPanel panel2;
    private JButton buttonAddRegion;
    private JButton button2;
    private JPanel panel3;
    private JCheckBox checkBoxZoomWhenNav;
    private JCheckBox checkBoxShowAllChrs;
    private CancelAction cancelAction;
    private AddRegionAction addRegionAction;
    private RemoveSelectedRegionsAction actionRemoveRegions;
    private ShowAllChromosomesCheckboxAction showAllChromosomesCheckboxAction;
    // JFormDesigner - End of variables declaration  //GEN-END:variables



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

    /**
     * Add a new RegionOfInterest for the current chromosome, with 0 start and end
     */
    private class AddRegionAction extends AbstractAction {
        private AddRegionAction() {
            // JFormDesigner - Action initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
            // Generated using JFormDesigner non-commercial license
            putValue(NAME, "Add Region");
            putValue(SHORT_DESCRIPTION, "Add a new region");
            // JFormDesigner - End of action initialization  //GEN-END:initComponents
        }

        public void actionPerformed(ActionEvent e) {
            String chr = FrameManager.getDefaultFrame().getChrName();
            if (chr == null || chr.isEmpty())
                JOptionPane.showMessageDialog(IGVMainFrame.getInstance(),
                        "No chromosome is specified. Can't create a region without a chromosome.",
                        "Error", JOptionPane.INFORMATION_MESSAGE);
            else if (chr.equalsIgnoreCase("All"))
                JOptionPane.showMessageDialog(IGVMainFrame.getInstance(),
                        "Regions cannot be created in the All Chromosomes view.",
                        "Error", JOptionPane.INFORMATION_MESSAGE);
            else
            {
                RegionOfInterest newRegion = new RegionOfInterest(FrameManager.getDefaultFrame().getChrName(),
                        0, 0, "");
                IGVMainFrame.getInstance().getSession().addRegionOfInterestWithNoListeners(newRegion);
            }
        }
    }

    private class RemoveSelectedRegionsAction extends AbstractAction {
        private RemoveSelectedRegionsAction() {
            // JFormDesigner - Action initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
            // Generated using JFormDesigner non-commercial license
            putValue(NAME, "Remove");
            putValue(SHORT_DESCRIPTION, "Remove all selected regions");
            // JFormDesigner - End of action initialization  //GEN-END:initComponents
        }

        public void actionPerformed(ActionEvent e) {
            int[] selectedRows = regionTable.getSelectedRows();
            if (selectedRows != null && selectedRows.length > 0)
            {
                List<RegionOfInterest> selectedRegions = new ArrayList<RegionOfInterest>();
                List<RegionOfInterest> regions = retrieveRegionsAsList();

                for (int selectedRowIndex : selectedRows)
                {
                    int selectedModelRow = regionTableRowSorter.convertRowIndexToModel(selectedRowIndex);
                    selectedRegions.add(regions.get(selectedModelRow));
                }
                regionTable.clearSelection();
                IGVMainFrame.getInstance().getSession().removeRegionsOfInterest(selectedRegions);

            }
            else
            {
                //todo dhmay -- I don't fully understand this call.  Clean this up.
                JOptionPane.showMessageDialog(IGVMainFrame.getInstance(), "No regions have been selected for removal.",
                        "Error", JOptionPane.INFORMATION_MESSAGE);
            }
        }
    }

    private class ShowAllChromosomesCheckboxAction extends AbstractAction {
        private ShowAllChromosomesCheckboxAction() {
            // JFormDesigner - Action initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
            // Generated using JFormDesigner non-commercial license
            putValue(NAME, "Show All Chrs");
            // JFormDesigner - End of action initialization  //GEN-END:initComponents
        }

        public void actionPerformed(ActionEvent e) {
            // TODO add your code here
            synchRegions();
        }
    }

    /**
     * Creates an appropriate popup for the row under the cursor, with Copy Sequence and Copy Details actions.
     * This class doesn't go back to the RegionOfInterest model -- it relies on the values stored in the
     * TableModel, since that's all we need.  It could easily grab the Region, though, like in RegionTableModelListener
     */
    private class RegionTablePopupHandler extends MouseAdapter {
        public void mousePressed(MouseEvent e) {
            if (SwingUtilities.isRightMouseButton(e)) {
                Point p = e.getPoint();
                //must convert row index from view to model, in case of sorting, filtering
                int row = regionTable.getRowSorter().convertRowIndexToModel(regionTable.rowAtPoint(p));
                int col = regionTable.columnAtPoint(p);

                if (row >= 0 && col >= 0) {
                    final String chr = (String) regionTableModel.getValueAt(row, TABLE_COLINDEX_CHR);
                    //displayed values are 1-based, so subract 1
                    final int start = (Integer) regionTableModel.getValueAt(row, TABLE_COLINDEX_START) - 1;
                    final int end = (Integer) regionTableModel.getValueAt(row, TABLE_COLINDEX_END) - 1;

                    final String desc = (String) regionTableModel.getValueAt(row, TABLE_COLINDEX_DESC);

                    JPopupMenu popupMenu = new JPopupMenu();
                    JMenuItem copySequenceItem = new JMenuItem("Copy Sequence");
                    copySequenceItem.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent e) {

                            String genomeId = GenomeManager.getInstance().getGenomeId();
                            byte[] seqBytes = SequenceManager.readSequence(genomeId, chr, start, end);
                            if (seqBytes == null) {
                                MessageUtils.showMessage("Sequence not available");
                            } else {
                                String sequence = new String(seqBytes);
                                copyTextToClipboard(sequence);
                            }
                        }
                    });
                    popupMenu.add(copySequenceItem);

                    JMenuItem copyDetailsItem = new JMenuItem("Copy Details");
                    copyDetailsItem.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent e) {
                            String details = chr + ":" + start + "-" + end;
                            if (desc != null && !desc.isEmpty())
                                details = details + ", " + desc;
                            copyTextToClipboard(details);
                        }
                    });
                    popupMenu.add(copyDetailsItem);
                    popupMenu.show(regionTable, p.x,p.y);
                }       
            }
        }

        /**
         * Copy a text string to the clipboard
         * @param text
         */
        private void copyTextToClipboard(String text)
        {
            StringSelection stringSelection = new StringSelection(text);
            Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
            clipboard.setContents(stringSelection, null);
        }
    }
}
