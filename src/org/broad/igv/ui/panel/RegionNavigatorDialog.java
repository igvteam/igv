/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Fred Hutchinson Cancer Research Center and Broad Institute
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


package org.broad.igv.ui.panel;

import com.google.common.eventbus.Subscribe;
import org.apache.log4j.Logger;
import org.broad.igv.feature.Range;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.lists.GeneList;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.event.ViewChange;
import org.broad.igv.util.StringUtils;

import javax.swing.*;
import javax.swing.border.LineBorder;
import javax.swing.event.*;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableColumnModel;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.util.List;

/**
 * @author Damon May
 *         <p/>
 *         This dialog displays a list of RegionOfInterest in a table and allows editing.
 *         Navigation to the start of a region is done by selecting the appropriate row.
 *         <p/>
 *         This dialog is not intended to be persistent.  To view one of these, create it.
 *         <p/>
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
     * Return the active RegionNavigatorDialog. null if none
     *
     * @return
     */
    public static RegionNavigatorDialog getOrCreateInstance(Frame owner) {
        if (activeInstance == null) {
            activeInstance = new RegionNavigatorDialog(owner);
        }
        return activeInstance;
    }

    /**
     * Return the active RegionNavigatorDialog. null if none
     *
     * @return
     */
    public static RegionNavigatorDialog getInstance() {
        return activeInstance;
    }

    /**
     * dispose the active instance and get rid of the pointer. Return whether or not there was an
     * active instance
     */
    public static boolean destroyInstance() {
        if (activeInstance == null)
            return false;
        activeInstance.dispose();
        activeInstance = null;
        return true;
    }

    private RegionNavigatorDialog(Frame owner) {
        super(owner);
        initComponents();
        postInit();
    }

    private RegionNavigatorDialog(Dialog owner) {
        super(owner);
        initComponents();
        postInit();
    }

    public void update(Observable observable, Object object) {
        synchRegions();
    }

    @Subscribe
    public void receiveChromosomeChanged(ViewChange.ChromosomeChangeResult e){
        synchRegions();
    }

    /**
     * Synchronize the regions ArrayList with the passed-in regionsCollection, and update UI
     */
    public void synchRegions() {
        //Indicate that we're synching regions, so that we don't respond to tableChanged events
        synchingRegions = true;
        List<RegionOfInterest> regions = retrieveRegionsAsList();
        regionTableModel = (DefaultTableModel) regionTable.getModel();
        while (regionTableModel.getRowCount() > 0){
            regionTableModel.removeRow(0);
        }
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

    private List<RegionOfInterest> retrieveRegionsAsList() {
        return new ArrayList<RegionOfInterest>(IGV.getInstance().getSession().getAllRegionsOfInterest());
    }

    /**
     * Populate the table with the loaded regions
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

        updateChromosomeDisplayed();

        synchRegions();

        ReferenceFrame defFrame = FrameManager.getDefaultFrame();
        defFrame.getEventBus().register(this);
        IGV.getInstance().getSession().getRegionsOfInterestObservable().addObserver(this);

        //resize window if small number of regions.  By default, tables are initialized with 20
        //rows, and that can look ungainly for empty windows or windows with a few rows.
        //This correction is rather hacky. Minimum size of 5 rows set.
        int newTableHeight = Math.min(regionTableModel.getRowCount() + 1, 5) * regionTable.getRowHeight();
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
        updateButtonsEnabled();
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
    public void updateChromosomeDisplayed() {
//        regionTable.updateUI();
//        showSearchedRegions();
        regionTableModel.fireTableDataChanged();
    }

    /**
     * Test whether we should display an entry
     *
     * @param regionChr
     * @param regionDesc
     * @return
     */
    protected boolean shouldIncludeRegion(String regionChr, String regionDesc) {
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
            if (firstRow > regions.size() - 1)
                return;
            //update all rows affected
            for (int i = firstRow; i <= Math.max(firstRow, Math.min(regionTable.getRowCount(), e.getLastRow())); i++)
                updateROIFromRegionTable(i);
        }
    }

    /**
     * Updates all ROIs with the values currently stored in the region table
     */
    public void updateROIsFromRegionTable() {
        for (int i = 0; i < regionTable.getRowSorter().getModelRowCount(); i++)
            updateROIFromRegionTable(i);
    }

    /**
     * Updates a single ROI with the values currently stored in the region table
     *
     * @param tableRow: the viewable index of the table row
     */
    public void updateROIFromRegionTable(int tableRow) {
        List<RegionOfInterest> regions = retrieveRegionsAsList();

        if (tableRow > regionTable.getRowCount() - 1)
            return;

        //must convert row index from view to model, in case of sorting, filtering
        int rowIdx = 0;

        try {
            rowIdx = regionTable.getRowSorter().convertRowIndexToModel(tableRow);
        } catch (ArrayIndexOutOfBoundsException x) {
            return;
        }

        RegionOfInterest region = regions.get(rowIdx);

        //dhmay changing 20110505: just update region values from all columns, instead of checking the event
        //to see which column is affected. This is in response to an intermittent bug.

        Object descObject = regionTableModel.getValueAt(rowIdx, TABLE_COLINDEX_DESC);
        if (descObject != null)
            region.setDescription(descObject.toString());

        //stored values are 0-based end-exclusive, viewed values are 1-based end-inclusive.  Check for negative number just in case
        int storeStartValue =
                Math.max(0, (Integer) regionTableModel.getValueAt(rowIdx, TABLE_COLINDEX_START) - 1);
        region.setStart(storeStartValue);

        int storeEndValue =
                Math.max(0, (Integer) regionTableModel.getValueAt(rowIdx, TABLE_COLINDEX_END));
        region.setEnd(storeEndValue);
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
                    for (int selectedRowIndex : selectedRows) {
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


//                    if (checkBoxZoomWhenNav.isSelected()) {
//                        // Option (1), zoom and center on group of selected regions, with an interval equal to
//                        // 20% of the length of the end regions on either side for context (dhmay reduced from 100%)
//                        int start = firstStartRegion.getStart() - (int) (0.2 * firstStartRegion.getLength());
//                        int end = lastEndRegion.getEnd() + (int) (0.2 * lastEndRegion.getLength());
//                        FrameManager.getDefaultFrame().jumpTo(selectedChrs.iterator().next(), start, end);
//                    } else {
//                        // Option (2), center on the FIRST selected region without changing resolution
//                        FrameManager.getDefaultFrame().centerOnLocation(firstStartRegion.getCenter());
//                    }

                }
            }
        }
    }


    /**
     * Return the selected regions in the table view.
     *
     * @param selectedRows
     * @return
     */
    private List<RegionOfInterest> getSelectedRegions(int[] selectedRows) {
        List<RegionOfInterest> selectedRegions = new ArrayList<RegionOfInterest>();
        List<RegionOfInterest> regions = retrieveRegionsAsList();

        for (int selectedRowIndex : selectedRows) {
            int selectedModelRow = regionTableRowSorter.convertRowIndexToModel(selectedRowIndex);
            selectedRegions.add(regions.get(selectedModelRow));
        }
        return selectedRegions;
    }

    private void regionTableMouseClicked(MouseEvent e) {
        //We update the state of associated buttons
        //whenever the user clicks a mouse button
        updateButtonsEnabled();
    }

    private void updateButtonsEnabled() {
        boolean enableMutStates = regionTable.getSelectedRowCount() >= 1;
        boolean enableZoomToRegion = regionTable.getSelectedRowCount() == 1;
        removeButton.setEnabled(enableMutStates);
        viewButton.setEnabled(enableMutStates);
        checkBoxZoomWhenNav.setEnabled(enableZoomToRegion);
    }

    private void thisWindowActivated(WindowEvent e) {
        synchRegions();
    }

    private void thisWindowDeactivated(WindowEvent e) {
        updateROIsFromRegionTable();
    }

    private void thisWindowClosed(WindowEvent e) {
        IGV.getInstance().getSession().getRegionsOfInterestObservable().deleteObserver(this);
        destroyInstance();
    }


    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        panel3 = new JPanel();
        checkBoxShowAllChrs = new JCheckBox();
        addButton = new JButton();
        removeButton = new JButton();
        scrollPane1 = new JScrollPane();
        regionTable = new JTable();
        panel1 = new JPanel();
        panel2 = new JPanel();
        viewButton = new JButton();
        checkBoxZoomWhenNav = new JCheckBox();
        panel4 = new JPanel();
        label1 = new JLabel();
        textFieldSearch = new JTextField();
        clearSearchButton = new JButton();
        cancelAction = new CancelAction();
        addAction = new AddRegionAction();
        actionRemoveRegions = new RemoveSelectedRegionsAction();
        showAllChromosomesCheckboxAction = new ShowAllChromosomesCheckboxAction();
        viewAction = new ViewSelectedAction();

        //======== this ========
        setTitle("Regions of Interest");
        setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        addWindowListener(new WindowAdapter() {
            @Override
            public void windowActivated(WindowEvent e) {
                thisWindowActivated(e);
            }
            @Override
            public void windowClosed(WindowEvent e) {
                thisWindowClosed(e);
            }
            @Override
            public void windowDeactivated(WindowEvent e) {
                thisWindowDeactivated(e);
            }
        });
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(null);
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(new BorderLayout());

                //======== panel3 ========
                {
                    panel3.setComponentOrientation(ComponentOrientation.LEFT_TO_RIGHT);
                    panel3.setBorder(LineBorder.createBlackLineBorder());
                    panel3.setLayout(new FlowLayout(FlowLayout.LEFT));

                    //---- checkBoxShowAllChrs ----
                    checkBoxShowAllChrs.setAction(showAllChromosomesCheckboxAction);
                    checkBoxShowAllChrs.setToolTipText("View regions from all chromosomes (otherwise, current chromosome only)");
                    checkBoxShowAllChrs.setSelected(true);
                    panel3.add(checkBoxShowAllChrs);

                    //---- addButton ----
                    addButton.setAction(addAction);
                    addButton.setText("Add");
                    addButton.setActionCommand("Add");
                    panel3.add(addButton);

                    //---- removeButton ----
                    removeButton.setAction(actionRemoveRegions);
                    removeButton.setText("Remove");
                    panel3.add(removeButton);
                }
                contentPanel.add(panel3, BorderLayout.NORTH);

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
                    regionTable.addMouseListener(new MouseAdapter() {
                        @Override
                        public void mouseClicked(MouseEvent e) {
                            regionTableMouseClicked(e);
                        }
                    });
                    scrollPane1.setViewportView(regionTable);
                }
                contentPanel.add(scrollPane1, BorderLayout.CENTER);

                //======== panel1 ========
                {
                    panel1.setLayout(new BoxLayout(panel1, BoxLayout.Y_AXIS));

                    //======== panel2 ========
                    {
                        panel2.setBorder(LineBorder.createBlackLineBorder());
                        panel2.setLayout(new FlowLayout(FlowLayout.LEFT));

                        //---- viewButton ----
                        viewButton.setText("View");
                        viewButton.setAction(viewAction);
                        viewButton.setActionCommand("View");
                        panel2.add(viewButton);

                        //---- checkBoxZoomWhenNav ----
                        checkBoxZoomWhenNav.setText("Zoom to Region");
                        checkBoxZoomWhenNav.setToolTipText("When navigating to a region, change zoom level?");
                        checkBoxZoomWhenNav.setSelected(true);
                        panel2.add(checkBoxZoomWhenNav);
                    }
                    panel1.add(panel2);

                    //======== panel4 ========
                    {
                        panel4.setBorder(LineBorder.createBlackLineBorder());
                        panel4.setLayout(new FlowLayout(FlowLayout.LEFT));

                        //---- label1 ----
                        label1.setText("Search");
                        panel4.add(label1);

                        //---- textFieldSearch ----
                        textFieldSearch.setToolTipText("Search for regions containing the specified description text.");
                        textFieldSearch.setPreferredSize(new Dimension(200, 28));
                        panel4.add(textFieldSearch);

                        //---- clearSearchButton ----
                        clearSearchButton.setAction(cancelAction);
                        clearSearchButton.setText("Clear Search");
                        panel4.add(clearSearchButton);
                    }
                    panel1.add(panel4);
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
    private JPanel panel3;
    private JCheckBox checkBoxShowAllChrs;
    private JButton addButton;
    private JButton removeButton;
    private JScrollPane scrollPane1;
    private JTable regionTable;
    private JPanel panel1;
    private JPanel panel2;
    private JButton viewButton;
    private JCheckBox checkBoxZoomWhenNav;
    private JPanel panel4;
    private JLabel label1;
    private JTextField textFieldSearch;
    private JButton clearSearchButton;
    private CancelAction cancelAction;
    private AddRegionAction addAction;
    private RemoveSelectedRegionsAction actionRemoveRegions;
    private ShowAllChromosomesCheckboxAction showAllChromosomesCheckboxAction;
    private ViewSelectedAction viewAction;
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
            putValue(NAME, "Add");
            putValue(SHORT_DESCRIPTION, "Add a new region");
            // JFormDesigner - End of action initialization  //GEN-END:initComponents
        }

        public void actionPerformed(ActionEvent e) {
            String chr = FrameManager.getDefaultFrame().getChrName();
            if (FrameManager.isGeneListMode()) {
                JOptionPane.showMessageDialog(IGV.getMainFrame(),
                        "Regions cannot be created in gene list or split-screen views.",
                        "Error", JOptionPane.INFORMATION_MESSAGE);

            } else if (chr == null || chr.isEmpty()) {
                JOptionPane.showMessageDialog(IGV.getMainFrame(),
                        "No chromosome is specified. Can't create a region without a chromosome.",
                        "Error", JOptionPane.INFORMATION_MESSAGE);
            } else if (chr.equalsIgnoreCase("All")) {
                JOptionPane.showMessageDialog(IGV.getMainFrame(),
                        "Regions cannot be created in the All Chromosomes view.",
                        "Error", JOptionPane.INFORMATION_MESSAGE);
            } else {
                Range r = FrameManager.getDefaultFrame().getCurrentRange();
                RegionOfInterest newRegion = new RegionOfInterest(r.getChr(), r.getStart(), r.getEnd(), "");
                IGV.getInstance().getSession().addRegionOfInterestWithNoListeners(newRegion);
            }
            updateButtonsEnabled();
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
            if (selectedRows != null && selectedRows.length > 0) {
                List<RegionOfInterest> selectedRegions = getSelectedRegions(selectedRows);
                IGV.getInstance().getSession().removeRegionsOfInterest(selectedRegions);
                synchRegions();
            } else {
                //todo dhmay -- I don't fully understand this call.  Clean this up.
                JOptionPane.showMessageDialog(IGV.getMainFrame(), "No regions have been selected for removal.",
                        "Error", JOptionPane.INFORMATION_MESSAGE);
            }
            updateButtonsEnabled();
        }
    }


    private class ViewSelectedAction extends AbstractAction {
        private ViewSelectedAction() {
            // JFormDesigner - Action initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
            // Generated using JFormDesigner non-commercial license
            putValue(NAME, "View");
            // JFormDesigner - End of action initialization  //GEN-END:initComponents
        }

        public void actionPerformed(ActionEvent e) {
            int[] selectedRows = regionTable.getSelectedRows();
            if (selectedRows != null && selectedRows.length > 0) {
                List<RegionOfInterest> selectedRegions = getSelectedRegions(selectedRows);

                // Create an "on-the-fly" gene list
                // TODO -- this is inefficient (converting regions -> strings then back again)
                List<String> loci = new ArrayList<String>(selectedRegions.size());
                if (checkBoxZoomWhenNav.isSelected() || selectedRegions.size() >= 2 || FrameManager.isGeneListMode()) {
                    for (RegionOfInterest roi : selectedRegions) {
                        loci.add(roi.getLocusString());
                    }
                } else {
                    //Need to preserve current zoom, iff checkbox not selected and only choosing 1
                    RegionOfInterest roi = selectedRegions.get(0);
                    Range range = FrameManager.getDefaultFrame().getCurrentRange();
                    int length = range.getLength();
                    int start = roi.getCenter() - length / 2;
                    int end = start + length;
                    //Shift so we don't go below 0
                    if(start < 0){
                        end += Math.abs(start);
                        start = 0;
                    }
                    loci.add(new RegionOfInterest(roi.getChr(), start, end, roi.getDescription()).getLocusString());
                }
                GeneList geneList = new GeneList("Regions of Interest", loci, false);
                IGV.getInstance().getSession().setCurrentGeneList(geneList);
                IGV.getInstance().resetFrames();

            }
            updateButtonsEnabled();
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
            updateButtonsEnabled();
        }
    }

    /**
     * Creates an appropriate popup for the row under the cursor, with Copy Sequence and Copy Details actions.
     * This class doesn't go back to the RegionOfInterest model -- it relies on the values stored in the
     * TableModel, since that's all we need.  It could easily grab the Region, though, like in RegionTableModelListener
     */
    private class RegionTablePopupHandler extends MouseAdapter {
        // Maximum length for "copy sequence" action
        private static final int MAX_SEQUENCE_LENGTH = 1000000;

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

                    JPopupMenu popupMenu = new IGVPopupMenu();
                    JMenuItem copySequenceItem = new JMenuItem("Copy Sequence");
                    copySequenceItem.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent e) {
                            int length = end - start;
                            if (length > MAX_SEQUENCE_LENGTH) {
                                JOptionPane.showMessageDialog(RegionNavigatorDialog.this, "Region is to large to copy sequence data.");
                            } else {
                                IGV.copySequenceToClipboard(GenomeManager.getInstance().getCurrentGenome(),
                                        chr, start, end);
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
                            StringUtils.copyTextToClipboard(details);
                        }
                    });
                    popupMenu.add(copyDetailsItem);
                    popupMenu.show(regionTable, p.x, p.y);
                }
            }
        }

    }

}
