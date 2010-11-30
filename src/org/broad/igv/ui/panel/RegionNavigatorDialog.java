/*
 * Created by JFormDesigner on Tue Nov 16 14:34:59 PST 2010
 */

package org.broad.igv.ui.panel;

import org.apache.log4j.Logger;
import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.ui.IGVMainFrame;

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

    //The active instance of RegionNavigatorDialog (only one at a time)
    public static RegionNavigatorDialog activeInstance;


    private DefaultTableModel regionTableModel;
    private List<RegionOfInterest> regions;

    private TableRowSorter<TableModel> regionTableRowSorter;

    public static RegionNavigatorDialog getActiveInstance() {
        return activeInstance;
    }

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
     * Synchronize the regions ArrayList with the passed-in regionsCollection, and update UI
     * @param regionsCollection
     */
    public void synchRegions(Collection<RegionOfInterest> regionsCollection)
    {
        regions = new ArrayList<RegionOfInterest>(regionsCollection != null ? regionsCollection :
                new ArrayList<RegionOfInterest>());
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




        regionTableModel.fireTableDataChanged();
    }

    /**
     * Populate the table with the loaded regions
     *
     * @param regionsCollection
     */
    private void postInit(Collection<RegionOfInterest> regionsCollection) {
        this.regions = new ArrayList<RegionOfInterest>(regionsCollection);
        regionTableModel = (DefaultTableModel) regionTable.getModel();


        regionTable.getSelectionModel().addListSelectionListener(new RegionTableSelectionListener());
        regionTableModel.addTableModelListener(new RegionTableModelListener());



        //custom row sorter required for displaying only a subset of rows
        regionTableRowSorter = new TableRowSorter<TableModel>(regionTableModel);
        regionTable.setRowSorter(regionTableRowSorter);
        regionTableRowSorter.setRowFilter(new RegionRowFilter());

        textFieldSearch.getDocument().addDocumentListener(new SearchFieldDocumentListener());

        checkBoxZoomWhenNav.setSelected(true);

        activeInstance = this;
        updateChromosomeDisplayed();

        synchRegions(regionsCollection);

        //resize window if small number of regions.  By default, tables are initialized with 20
        //rows, and that can look ungainly for empty windows or windows with a few rows.
        //This correction is rather hacky. Minimum size of 5 rows set.
        int newTableHeight = Math.min(regions.size()+1,5) * regionTable.getRowHeight();
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
     * A row filter that shows only rows that contain filterString, case-insensitive
     */
    private class RegionRowFilter extends RowFilter<TableModel, Object> {

        public RegionRowFilter() {
            super();
        }

        public boolean include(RowFilter.Entry entry) {
            //if table is empty, a non-region event is fed here.  Test for it and don't display
            if (entry.getValue(TABLE_COLINDEX_CHR) == null)
                return false;
            
            String filterStringLowercase = null;
            if (textFieldSearch.getText() != null)
                filterStringLowercase = textFieldSearch.getText().toLowerCase();

            String regionDesc = (String) entry.getValue(TABLE_COLINDEX_DESC);
            String chr = FrameManager.getDefaultFrame().getChrName();
            //show only regions in the current chromosome
            if (chr != null && !chr.isEmpty() && !entry.getValue(TABLE_COLINDEX_CHR).equals(chr))
                return false;
            //show only regions matching the search string (if specified)
            if ((filterStringLowercase != null && !filterStringLowercase.isEmpty() &&
                (regionDesc == null || !regionDesc.toLowerCase().contains(filterStringLowercase))))
                return false;

            return true;
        }
    }

    /**
     * Listen for updates to the cells, save changes to the Regions
     */
    private class RegionTableModelListener implements TableModelListener {
        public void tableChanged(TableModelEvent e) {
            int rowIdx = e.getFirstRow();
            //range checking because this method gets called after a clear event
            if (rowIdx > regions.size()-1)
                return;
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
                int[] selectedRows = regionTable.getSelectedRows();
                if (selectedRows != null && selectedRows.length > 0
                        && regions.size() >= selectedRows.length)  //dhmay: this is hacky. Bad things can happen with clear regions
                {
                    RegionOfInterest firstStartRegion = null;
                    RegionOfInterest lastEndRegion = null;

                    //Figure out which region has the first start and which has the last end
                    for (int selectedRowIndex : selectedRows)
                    {
                        int selectedModelRow = regionTableRowSorter.convertRowIndexToModel(selectedRowIndex);
                        RegionOfInterest region = regions.get(selectedModelRow);
                        if (firstStartRegion == null || region.getStart() < firstStartRegion.getStart())
                            firstStartRegion = region;
                        if (lastEndRegion == null || region.getEnd() > lastEndRegion.getEnd())
                            lastEndRegion = region;
                    }

                    if (checkBoxZoomWhenNav.isSelected())
                    {
                        // Option (1), zoom and center on group of selected regions, with an interval equal to
                        // 20% of the length of the end regions on either side for context (dhmay reduced from 100%)
                        int start = firstStartRegion.getStart() - (int) (0.2 * firstStartRegion.getLength());
                        int end = lastEndRegion.getEnd() + (int) (0.2 * lastEndRegion.getLength());
                        FrameManager.getDefaultFrame().jumpTo(FrameManager.getDefaultFrame().getChrName(), start, end);
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
        // Generated using JFormDesigner Evaluation license - Damon May
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        scrollPane1 = new JScrollPane();
        regionTable = new JTable();
        panel1 = new JPanel();
        textFieldSearch = new JTextField();
        label1 = new JLabel();
        button1 = new JButton();
        panel2 = new JPanel();
        checkBoxZoomWhenNav = new JCheckBox();
        buttonAddRegion = new JButton();
        button2 = new JButton();
        cancelAction = new CancelAction();
        addRegionAction = new AddRegionAction();
        actionRemoveRegions = new RemoveSelectedRegionsAction();

        //======== this ========
        setTitle("Regions of Interest");
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(null);

            // JFormDesigner evaluation mark
            dialogPane.setBorder(new javax.swing.border.CompoundBorder(
                new javax.swing.border.TitledBorder(new javax.swing.border.EmptyBorder(0, 0, 0, 0),
                    "JFormDesigner Evaluation", javax.swing.border.TitledBorder.CENTER,
                    javax.swing.border.TitledBorder.BOTTOM, new java.awt.Font("Dialog", java.awt.Font.BOLD, 12),
                    java.awt.Color.red), dialogPane.getBorder())); dialogPane.addPropertyChangeListener(new java.beans.PropertyChangeListener(){public void propertyChange(java.beans.PropertyChangeEvent e){if("border".equals(e.getPropertyName()))throw new RuntimeException();}});

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

                        //---- checkBoxZoomWhenNav ----
                        checkBoxZoomWhenNav.setText("Zoom when Navigating");
                        panel2.add(checkBoxZoomWhenNav, BorderLayout.EAST);

                        //---- buttonAddRegion ----
                        buttonAddRegion.setAction(addRegionAction);
                        panel2.add(buttonAddRegion, BorderLayout.WEST);

                        //---- button2 ----
                        button2.setAction(actionRemoveRegions);
                        button2.setText("Remove Selected");
                        panel2.add(button2, BorderLayout.CENTER);
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
    // Generated using JFormDesigner Evaluation license - Damon May
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JScrollPane scrollPane1;
    private JTable regionTable;
    private JPanel panel1;
    private JTextField textFieldSearch;
    private JLabel label1;
    private JButton button1;
    private JPanel panel2;
    private JCheckBox checkBoxZoomWhenNav;
    private JButton buttonAddRegion;
    private JButton button2;
    private CancelAction cancelAction;
    private AddRegionAction addRegionAction;
    private RemoveSelectedRegionsAction actionRemoveRegions;
    // JFormDesigner - End of variables declaration  //GEN-END:variables



    private class CancelAction extends AbstractAction {
        private CancelAction() {
            // JFormDesigner - Action initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
            // Generated using JFormDesigner Evaluation license - Damon May
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
            // Generated using JFormDesigner Evaluation license - Damon May
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
            // Generated using JFormDesigner Evaluation license - Damon May
            putValue(NAME, "Remove");
            putValue(SHORT_DESCRIPTION, "Remove all selected regions");
            // JFormDesigner - End of action initialization  //GEN-END:initComponents
        }

        public void actionPerformed(ActionEvent e) {
            int[] selectedRows = regionTable.getSelectedRows();
            if (selectedRows != null && selectedRows.length > 0)
            {
                List<RegionOfInterest> selectedRegions = new ArrayList<RegionOfInterest>();

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
}
