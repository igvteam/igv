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
 * Created by JFormDesigner on Tue Mar 06 14:09:15 EST 2012
 */

package org.broad.igv.cbio;

import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.PreferenceManager;
import org.broad.igv.lists.GeneList;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.ui.util.*;
import org.broad.igv.util.BrowserLauncher;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.StringUtils;
import org.w3c.dom.Node;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.border.EtchedBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;
import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.List;

/**
 * Dialog for letting the user filter a GeneNetwork from
 * cBio.
 *
 * @author Jacob Silterra
 */
public class FilterGeneNetworkUI extends JDialog {

    private Logger log = Logger.getLogger(FilterGeneNetworkUI.class);

    private GeneList seedGeneList;
    private List<AttributeFilter> filterRows = new ArrayList<AttributeFilter>(1);
    GeneNetwork network = null;

    private GraphListModel listModel;
    private Map<String, JTextField> thresholdsMap = new HashMap<String, JTextField>(5);


    private static List<String> columnNames;
    private static Map<Integer, String> columnNumToKeyMap;

    /**
     * For change events, we keep track of where we started from
     */
    private Component lastSelectedTab;

    static {
        String[] firstLabels = {"Gene label", "Interactions"};
        columnNumToKeyMap = new HashMap<Integer, String>(GeneNetwork.attributeMap.size());
        columnNames = new ArrayList<String>(firstLabels.length + GeneNetwork.attributeMap.size());
        for (String label : firstLabels) {
            columnNames.add(label);
        }
        int ind = columnNames.size();
        for (String key : GeneNetwork.attributeMap.keySet()) {
            columnNumToKeyMap.put(ind, key);
            ind++;
            columnNames.add(AttributeFilter.keyToLabel(key));

        }

    }

    public FilterGeneNetworkUI(Frame owner, GeneList seedGeneList) {
        super(owner);
        this.seedGeneList = seedGeneList;
    }

    @Override
    public void setVisible(boolean visible) {
        if (visible && network == null) {
            loadcBioData();
        } else {
            super.setVisible(visible);
        }
    }

    private void loadcBioData() {
        network = null;
        final List<String> geneLoci = seedGeneList.getLoci();
        final IndefiniteProgressMonitor indefMonitor = new IndefiniteProgressMonitor();
        final ProgressBar.ProgressDialog progressDialog = ProgressBar.showProgressDialog((Frame) getOwner(), "Loading cBio data...", indefMonitor, true);
        progressDialog.getProgressBar().setIndeterminate(true);
        indefMonitor.start();

        //Since we load the data asynchronously, we run this when finished
        final Runnable updateUI = new Runnable() {
            @Override
            public void run() {
                boolean needSetup = !FilterGeneNetworkUI.this.isVisible();
                if (needSetup) {
                    initComponents();
                }

                setGeneNetworkNeedsUpdate(false);
                tabbedPane.setSelectedComponent(filterPane);
                initComponentData();

                if (needSetup) {
                    FilterGeneNetworkUI.super.setVisible(true);
                }
            }
        };

        final Runnable runnable = new Runnable() {
            @Override
            public void run() {

                WaitCursorManager.CursorToken token = null;

                try {
                    token = WaitCursorManager.showWaitCursor();
                    network = GeneNetwork.getFromCBIO(geneLoci);
                    if (network.vertexSet().size() == 0) {
                        MessageUtils.showMessage("No results found for " + StringUtils.join(geneLoci, ", "));
                    } else {
                        network.annotateAll(IGV.getInstance().getAllTracks());
                        UIUtilities.invokeOnEventThread(updateUI);
                    }
                } catch (Throwable e) {
                    e.printStackTrace();
                    log.error(e.getMessage());
                    MessageUtils.showMessage("Error loading data: " + e.getMessage());
                } finally {
                    WaitCursorManager.removeWaitCursor(token);

                    if (progressDialog != null) {
                        progressDialog.setVisible(false);
                        indefMonitor.stop();
                    }
                }
            }
        };

        LongRunningTask.submit(runnable);
    }


    private void initThresholdsMap() {
        thresholdsMap.put(PreferenceManager.CBIO_MUTATION_THRESHOLD, mutInput);
        thresholdsMap.put(PreferenceManager.CBIO_AMPLIFICATION_THRESHOLD, ampInput);
        thresholdsMap.put(PreferenceManager.CBIO_DELETION_THRESHOLD, delInput);
        thresholdsMap.put(PreferenceManager.CBIO_EXPRESSION_UP_THRESHOLD, expUpInput);
        thresholdsMap.put(PreferenceManager.CBIO_EXPRESSION_DOWN_THRESHOLD, expDownInput);
    }

    /**
     * Load data into visual components. This should be called
     * AFTER loading network.
     */
    private void initComponentData() {
        refreshSeedGenesTextArea();

        if (this.filterRows.size() == 0) {
            add();
        }

        initThresholdsMap();
        loadThresholds();

        listModel = new GraphListModel();
        geneTable.setModel(listModel);
        initRenderers();

        applySoftFilters();
    }

    private void modifyForSelection(JComponent component, boolean isSelected) {
        if (isSelected) {
            component.setForeground(geneTable.getSelectionForeground());
            component.setBackground(geneTable.getSelectionBackground());
            component.setOpaque(isSelected);
        }
    }

    private void initRenderers() {
        //Bold seed genes and format numbers nicely
        TableCellRenderer stringRenderer = new TableCellRenderer() {
            @Override
            public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
                JLabel comp = new JLabel();
                if (value != null) {
                    boolean isSeedGene = seedGeneList.getLoci().contains(value);
                    comp.setText(String.valueOf(value));
                    if (isSeedGene) {
                        comp.setFont(comp.getFont().deriveFont(comp.getFont().getStyle() | Font.BOLD));
                    }
                }
                modifyForSelection(comp, isSelected);
                return comp;
            }
        };

        TableCellRenderer intRenderer = new TableCellRenderer() {
            @Override
            public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
                JLabel comp = new JLabel("" + value);
                modifyForSelection(comp, isSelected);
                return comp;
            }
        };


        TableCellRenderer doubleRenderer = new TableCellRenderer() {
            @Override
            public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {

                JLabel comp = new JLabel();
                if (value != null) {
                    String sVal;
                    Double dPerc = (Double) value;
                    if (dPerc == 0.0d) {
                        sVal = "0.0";
                    } else {
                        //If small, show in exponential format
                        //Otherwise just show 1 decimal place
                        String fmt = "%2.1f";
                        if (dPerc < 0.1d) {
                            fmt = "%2.1e";
                        }
                        sVal = String.format(fmt, dPerc);
                    }
                    comp.setText(sVal);
                }
                modifyForSelection(comp, isSelected);
                return comp;
            }
        };

        geneTable.setDefaultRenderer(String.class, stringRenderer);
        geneTable.setDefaultRenderer(Integer.class, intRenderer);
        geneTable.setDefaultRenderer(Double.class, doubleRenderer);
    }

    /**
     * Update the displayed text area to show the genes in {@code seedGeneList}
     * See {@linkplain #getNewSeedGeneList}
     */
    private void refreshSeedGenesTextArea() {
        String seedGenesString = StringUtils.join(seedGeneList.getLoci(), "\n");
        seedGenesText.setText(seedGenesString);
    }

    /**
     * Get the text content of {@code seedGenestext} and parse it into
     * a GeneList.
     * See {@linkplain #refreshSeedGenesTextArea()}
     *
     * @return new seedGeneList, which may be the same as the old
     */
    private GeneList getNewSeedGeneList() {
        String[] genes = seedGenesText.getText().toUpperCase().split("[\\r\\n]{1,2}");
        boolean updated = !Arrays.equals(seedGeneList.getLoci().toArray(), genes);
        if (updated) {
            return new GeneList(null, Arrays.asList(genes));
        } else {
            return seedGeneList;
        }
    }

    private boolean checkNewSeedGeneList() {
        return getNewSeedGeneList() != seedGeneList;
    }

    private void remove(AttributeFilter row) {
        contentPane.remove(row.getPanel());
        filterRows.remove(row);
        int numRows = filterRows.size();
        filterRows.get(numRows - 1).setIsLast(true);
        filterRows.get(0).setShowDel(numRows >= 2);

        validate();

        adjustWindowHeight(-row.getPanel().getHeight());

        applySoftFilters();
    }


    private void add() {
        final AttributeFilter row = new AttributeFilter();

        row.getDelRow().addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                remove(row);
            }
        });

        row.getAddRow().addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                add();
            }
        });
        contentPane.add(row.getPanel());

        //We want to refresh the
        RefreshListener listener = new RefreshListener();

        row.getAttrName().addActionListener(listener);
        for (JTextField text : new JTextField[]{row.minVal, row.maxVal}) {
            text.addActionListener(listener);
            text.addFocusListener(listener);
        }


        //Set the status of being last
        if (filterRows.size() >= 1) {
            filterRows.get(filterRows.size() - 1).setIsLast(false);
        }
        filterRows.add(row);

        int numRows = filterRows.size();
        filterRows.get(numRows - 1).setIsLast(true);
        filterRows.get(0).setShowDel(numRows >= 2);

        validate();

        adjustWindowHeight(row.getPanel().getHeight());
    }

    /**
     * When adding/removing rows, we change the height of some components
     * This method may have no effect
     *
     * @param heightDelta Number of pixels to increase (negative values allowed)
     *                    the tab pane and dialog itself
     */
    private void adjustWindowHeight(int heightDelta) {

        int newWinHeight = getHeight() + heightDelta;
        if (filterRows.size() >= 4 && (newWinHeight < Toolkit.getDefaultToolkit().getScreenSize().getHeight() - 100)) {
            Dimension newSize = tabbedPane.getPreferredSize();
            newSize.setSize(newSize.getWidth(), newSize.getHeight() + heightDelta);
            tabbedPane.setPreferredSize(newSize);

            newSize = getSize();
            newSize.setSize(newSize.getWidth(), newWinHeight);
            setSize(newSize);

            validate();
        }

    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        setVisible(false);
    }

    private void applySoftFilters() {
        network.reset();

        if (showSeedOnly.isSelected()) {
            //Redundant, as filterGenes already preserves query genes
            //But we want to be explicit
            network.filterGenes(GeneNetwork.inQuery);
        } else {
            //TODO This is only AND, should also include OR
            for (AttributeFilter filter : this.filterRows) {
                String filt_el = (String) filter.getAttrName().getSelectedItem();
                if (GeneNetwork.attributeMap.containsKey(filt_el) || GeneNetwork.PERCENT_ALTERED.equals(filt_el)) {
                    float min = Float.parseFloat(filter.minVal.getText());
                    float max = Float.parseFloat(filter.maxVal.getText());
                    network.filterGenesRange(filt_el, min / 100, max / 100);
                }
            }
            if (!keepIsolated.isSelected()) {
                network.pruneGraph();
            }
        }

        totNumGenes.setText("Total Genes: " + network.geneVertexes().size());

        this.listModel.markDirty();
    }

    private void showNetwork() {

        Runnable runnable = new Runnable() {
            @Override
            public void run() {
                try {
                    String url = network.outputForcBioView();
                    url = "file://" + url;
                    BrowserLauncher.openURL(url);
                } catch (IOException err) {
                    log.error(err);
                    MessageUtils.showMessage("Error opening network for viewing. " + err.getMessage());
                }
            }
        };
        LongRunningTask.submit(runnable);
    }

    /**
     * Sorting of the rows is done on the view, not the underlying model,
     * so we must convert. See {@link TableRowSorter}
     *
     * @return int[] of model indices which are selected
     */
    private int[] getModelIndices() {
        int[] selection = geneTable.getSelectedRows();
        for (int i = 0; i < selection.length; i++) {
            selection[i] = geneTable.convertRowIndexToModel(selection[i]);
        }
        return selection;
    }

    private void okButtonActionPerformed(ActionEvent e) {
        //setVisible(false);
        showNetwork();
    }

    private void addRowActionPerformed(ActionEvent e) {
        add();
    }

    /**
     * Refresh the view based on filter input.
     */
    private void refreshFilters() {
        this.applySoftFilters();
        listModel.markDirty();
        this.validate();
    }

    private void refFilterActionPerformed(ActionEvent e) {
        refreshFilters();
    }

    private boolean saveThresholds() {
        try {
            for (Map.Entry<String, JTextField> entry : thresholdsMap.entrySet()) {
                float fval = Float.parseFloat(entry.getValue().getText());
                String sval = "" + fval;
                if (entry.getKey() == PreferenceManager.CBIO_MUTATION_THRESHOLD) {
                    int ival = Integer.parseInt(entry.getValue().getText());
                    sval = "" + ival;
                }
                PreferenceManager.getInstance().put(entry.getKey(), sval);
            }
        } catch (NumberFormatException e) {
            MessageUtils.showMessage("Invalid input. Mutation count must be integers, others must be numeric. " + e.getMessage());
            return false;
        }
        return true;
    }

    private void loadThresholds() {
        for (Map.Entry<String, JTextField> entry : thresholdsMap.entrySet()) {
            String value = PreferenceManager.getInstance().get(entry.getKey());
            entry.getValue().setText(value);
        }
    }

    private void tabbedPaneStateChanged(ChangeEvent e) {
        if (tabbedPane == null || thresholdsPane == null) {
            //Component not built yet
            return;
        }

        if (lastSelectedTab != null) {

            //If the user enters an invalid threshold, we don't let them switch away
            if (lastSelectedTab == seedGenesPane && !saveThresholds()) {
                tabbedPane.setSelectedComponent(thresholdsPane);
            }

            //Reload network when we change away from seed genes pane
            //No-op if unchanged
//            if (lastSelectedTab == seedGenesPane) {
//                updateNetwork();
//            }
        }
        lastSelectedTab = tabbedPane.getSelectedComponent();
    }

    private void saveButtonActionPerformed(ActionEvent e) {
        File outPath = FileDialogUtils.chooseFile("Save table to...", DirectoryManager.getUserDirectory(), FileDialogUtils.SAVE);
        if (outPath != null) {
            try {
                saveTable(outPath);
            } catch (FileNotFoundException exc) {
                MessageUtils.showMessage(exc.getMessage());
            }
        }
    }

    private void resetToDefaultsButtonActionPerformed(ActionEvent e) {
        for (Map.Entry<String, JTextField> entry : thresholdsMap.entrySet()) {
            String value = PreferenceManager.getInstance().getDefaultValue(entry.getKey());
            entry.getValue().setText(value);
        }
        saveThresholds();
    }

    private void cancel2ActionPerformed(ActionEvent e) {
        cancelButton.doClick();
    }

    private void retrieveNetworkButtonActionPerformed(ActionEvent e) {
        updateNetwork();
    }

    private void updateNetwork() {
        GeneList newSeedGeneList = getNewSeedGeneList();
        boolean updated = newSeedGeneList != seedGeneList;
        if (updated) {
            seedGeneList = newSeedGeneList;
            loadcBioData();
        }
    }

    private void resetSeedGeneTextAreaButtonActionPerformed(ActionEvent e) {
        refreshSeedGenesTextArea();
        setGeneNetworkNeedsUpdate(false);
    }

    private void seedGenesTextKeyReleased(KeyEvent e) {
        boolean updated = checkNewSeedGeneList();
        setGeneNetworkNeedsUpdate(updated);
    }

    private void setGeneNetworkNeedsUpdate(boolean updated) {
        boolean tabsEnabled = !updated;
        int filterIndex = tabbedPane.indexOfComponent(filterPane);
        int thresholdsIndex = tabbedPane.indexOfComponent(thresholdsPane);
        int[] indexes = {filterIndex, thresholdsIndex};
        Color color = tabsEnabled ? Color.black : Color.gray;
        for (int index : indexes) {
            tabbedPane.setEnabledAt(index, tabsEnabled);
            tabbedPane.setForegroundAt(index, color);
        }

        retrieveNetworkButton.setEnabled(updated);
    }

    private void showSeedOnlyActionPerformed(ActionEvent e) {
        refreshFilters();
    }


    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        tabbedPane = new JTabbedPane();
        seedGenesPane = new JPanel();
        vSpacer1 = new JPanel(null);
        label5 = new JLabel();
        scrollPane2 = new JScrollPane();
        seedGenesText = new JTextArea();
        textPane1 = new JTextPane();
        seedButtonBar = new JPanel();
        retrieveNetworkButton = new JButton();
        resetSeedGeneTextAreaButton = new JButton();
        cancel2 = new JButton();
        filterPane = new JPanel();
        panel1 = new JPanel();
        addRow = new JButton();
        contentPane = new JPanel();
        scrollPane1 = new JScrollPane();
        geneTable = new JTable();
        buttonBar = new JPanel();
        totNumGenes = new JLabel();
        showSeedOnly = new JCheckBox();
        refFilter = new JButton();
        keepIsolated = new JCheckBox();
        okButton = new JButton();
        cancelButton = new JButton();
        saveButton = new JButton();
        helpButton = new JButton();
        thresholdsPane = new JPanel();
        contentPanel = new JPanel();
        label2 = new JLabel();
        label3 = new JLabel();
        delInput = new JTextField();
        label4 = new JLabel();
        expUpInput = new JTextField();
        label7 = new JLabel();
        expDownInput = new JTextField();
        ampInput = new JTextField();
        label1 = new JLabel();
        mutInput = new JTextField();
        label6 = new JLabel();
        label8 = new JLabel();
        separator1 = new JSeparator();
        separator3 = new JSeparator();
        separator2 = new JSeparator();
        panel2 = new JPanel();
        textArea1 = new JTextArea();
        resetToDefaultsButton = new JButton();

        //======== this ========
        setMinimumSize(new Dimension(600, 22));
        setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        setModal(true);
        setModalityType(Dialog.ModalityType.DOCUMENT_MODAL);
        Container contentPane2 = getContentPane();
        contentPane2.setLayout(new BorderLayout());

        //======== tabbedPane ========
        {
            tabbedPane.setPreferredSize(new Dimension(571, 400));
            tabbedPane.setMinimumSize(new Dimension(571, 346));
            tabbedPane.addChangeListener(new ChangeListener() {
                @Override
                public void stateChanged(ChangeEvent e) {
                    tabbedPaneStateChanged(e);
                }
            });

            //======== seedGenesPane ========
            {
                seedGenesPane.setLayout(new GridBagLayout());
                ((GridBagLayout) seedGenesPane.getLayout()).rowHeights = new int[]{0, 0, 0, 0, 0, 0};
                ((GridBagLayout) seedGenesPane.getLayout()).rowWeights = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 1.0E-4};
                seedGenesPane.add(vSpacer1, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 5, 0), 0, 0));

                //---- label5 ----
                label5.setText("Seed Genes (one per line):");
                label5.setLabelFor(seedGenesText);
                seedGenesPane.add(label5, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 5, 0), 0, 0));

                //======== scrollPane2 ========
                {
                    scrollPane2.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);

                    //---- seedGenesText ----
                    seedGenesText.setRows(12);
                    seedGenesText.setToolTipText("cBio will be queried to find what genes interact with these genes");
                    seedGenesText.setDragEnabled(false);
                    seedGenesText.addKeyListener(new KeyAdapter() {
                        @Override
                        public void keyReleased(KeyEvent e) {
                            seedGenesTextKeyReleased(e);
                        }
                    });
                    scrollPane2.setViewportView(seedGenesText);
                }
                seedGenesPane.add(scrollPane2, new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 5, 0), 0, 0));

                //---- textPane1 ----
                textPane1.setBorder(null);
                textPane1.setEditable(false);
                textPane1.setText("IGV will query cBio to find genes which interact with the seed genes entered here.");
                textPane1.setBackground(UIManager.getColor("Button.background"));
                seedGenesPane.add(textPane1, new GridBagConstraints(0, 3, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 5, 0), 0, 0));

                //======== seedButtonBar ========
                {
                    seedButtonBar.setLayout(new GridBagLayout());
                    ((GridBagLayout) seedButtonBar.getLayout()).columnWidths = new int[]{0, 0, 0, 0};
                    ((GridBagLayout) seedButtonBar.getLayout()).rowHeights = new int[]{0, 0};
                    ((GridBagLayout) seedButtonBar.getLayout()).columnWeights = new double[]{0.0, 0.0, 0.0, 1.0E-4};
                    ((GridBagLayout) seedButtonBar.getLayout()).rowWeights = new double[]{0.0, 1.0E-4};

                    //---- retrieveNetworkButton ----
                    retrieveNetworkButton.setText("Retrieve Network");
                    retrieveNetworkButton.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            retrieveNetworkButtonActionPerformed(e);
                        }
                    });
                    seedButtonBar.add(retrieveNetworkButton, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0,
                            GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                            new Insets(0, 0, 0, 5), 0, 0));

                    //---- resetSeedGeneTextAreaButton ----
                    resetSeedGeneTextAreaButton.setText("Reset to Original");
                    resetSeedGeneTextAreaButton.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            resetSeedGeneTextAreaButtonActionPerformed(e);
                        }
                    });
                    seedButtonBar.add(resetSeedGeneTextAreaButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
                            GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                            new Insets(0, 0, 0, 5), 0, 0));

                    //---- cancel2 ----
                    cancel2.setText("Cancel");
                    cancel2.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            cancel2ActionPerformed(e);
                        }
                    });
                    seedButtonBar.add(cancel2, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0,
                            GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                            new Insets(0, 0, 0, 0), 0, 0));
                }
                seedGenesPane.add(seedButtonBar, new GridBagConstraints(0, 4, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 0, 0), 0, 0));
            }
            tabbedPane.addTab("Seed Genes", seedGenesPane);


            //======== filterPane ========
            {
                filterPane.setBorder(new EmptyBorder(12, 12, 12, 12));
                filterPane.setMinimumSize(new Dimension(443, 300));
                filterPane.setLayout(new GridBagLayout());
                ((GridBagLayout) filterPane.getLayout()).columnWidths = new int[]{0, 0};
                ((GridBagLayout) filterPane.getLayout()).rowHeights = new int[]{0, 0, 0, 0, 0, 0};
                ((GridBagLayout) filterPane.getLayout()).columnWeights = new double[]{1.0, 1.0E-4};
                ((GridBagLayout) filterPane.getLayout()).rowWeights = new double[]{0.0, 0.0, 1.0, 0.0, 0.0, 1.0E-4};

                //======== panel1 ========
                {
                    panel1.setLayout(new GridBagLayout());
                    ((GridBagLayout) panel1.getLayout()).columnWidths = new int[]{0, 0, 0};
                    ((GridBagLayout) panel1.getLayout()).rowHeights = new int[]{0, 0};
                    ((GridBagLayout) panel1.getLayout()).columnWeights = new double[]{0.0, 0.0, 1.0E-4};
                    ((GridBagLayout) panel1.getLayout()).rowWeights = new double[]{0.0, 1.0E-4};

                    //---- addRow ----
                    addRow.setText("Add Filter");
                    addRow.setMaximumSize(new Dimension(200, 28));
                    addRow.setMinimumSize(new Dimension(100, 28));
                    addRow.setPreferredSize(new Dimension(150, 28));
                    addRow.setVisible(false);
                    addRow.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            addRowActionPerformed(e);
                        }
                    });
                    panel1.add(addRow, new GridBagConstraints(0, 0, 2, 1, 0.0, 0.0,
                            GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                            new Insets(0, 0, 0, 0), 0, 0));
                }
                filterPane.add(panel1, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 0, 0), 0, 0));

                //======== contentPane ========
                {
                    contentPane.setLayout(new BoxLayout(contentPane, BoxLayout.Y_AXIS));
                }
                filterPane.add(contentPane, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 0, 0), 0, 0));

                //======== scrollPane1 ========
                {

                    //---- geneTable ----
                    geneTable.setAutoCreateRowSorter(true);
                    scrollPane1.setViewportView(geneTable);
                }
                filterPane.add(scrollPane1, new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 0, 0), 0, 0));

                //======== buttonBar ========
                {
                    buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                    buttonBar.setMaximumSize(new Dimension(2147483647, 137));
                    buttonBar.setPreferredSize(new Dimension(421, 100));
                    buttonBar.setMinimumSize(new Dimension(421, 80));
                    buttonBar.setLayout(new GridBagLayout());
                    ((GridBagLayout) buttonBar.getLayout()).columnWidths = new int[]{0, 85, 85, 80};
                    ((GridBagLayout) buttonBar.getLayout()).columnWeights = new double[]{1.0, 0.0, 0.0, 0.0};

                    //---- totNumGenes ----
                    totNumGenes.setText("Total Genes: #");
                    buttonBar.add(totNumGenes, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0,
                            GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                            new Insets(0, 0, 5, 5), 0, 0));

                    //---- showSeedOnly ----
                    showSeedOnly.setText("Seed Genes Only");
                    showSeedOnly.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            showSeedOnlyActionPerformed(e);
                        }
                    });
                    buttonBar.add(showSeedOnly, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
                            GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                            new Insets(0, 0, 5, 5), 0, 0));

                    //---- refFilter ----
                    refFilter.setText("Refresh Filter");
                    refFilter.setVisible(false);
                    refFilter.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            refFilterActionPerformed(e);
                        }
                    });
                    buttonBar.add(refFilter, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0,
                            GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                            new Insets(0, 0, 5, 5), 0, 0));

                    //---- keepIsolated ----
                    keepIsolated.setText("Keep Isolated Genes");
                    keepIsolated.setVisible(false);
                    buttonBar.add(keepIsolated, new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0,
                            GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                            new Insets(0, 0, 5, 0), 0, 0));

                    //---- okButton ----
                    okButton.setText("View Network");
                    okButton.setToolTipText("Display the network in a web browser");
                    okButton.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            okButtonActionPerformed(e);
                        }
                    });
                    buttonBar.add(okButton, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0,
                            GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                            new Insets(0, 0, 0, 5), 0, 0));

                    //---- cancelButton ----
                    cancelButton.setText("Cancel");
                    cancelButton.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            cancelButtonActionPerformed(e);
                        }
                    });
                    buttonBar.add(cancelButton, new GridBagConstraints(3, 1, 1, 1, 0.0, 0.0,
                            GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                            new Insets(0, 0, 0, 0), 0, 0));

                    //---- saveButton ----
                    saveButton.setText("Save Table");
                    saveButton.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            saveButtonActionPerformed(e);
                        }
                    });
                    buttonBar.add(saveButton, new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0,
                            GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                            new Insets(0, 0, 0, 5), 0, 0));

                    //---- helpButton ----
                    helpButton.setText("Help");
                    helpButton.setVisible(false);
                    buttonBar.add(helpButton, new GridBagConstraints(2, 1, 1, 1, 0.0, 0.0,
                            GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                            new Insets(0, 0, 0, 5), 0, 0));
                }
                filterPane.add(buttonBar, new GridBagConstraints(0, 3, 1, 1, 0.0, 0.0,
                        GridBagConstraints.NORTH, GridBagConstraints.HORIZONTAL,
                        new Insets(0, 0, 0, 0), 0, 0));
            }
            tabbedPane.addTab("Filter", filterPane);


            //======== thresholdsPane ========
            {
                thresholdsPane.setPreferredSize(new Dimension(550, 196));
                thresholdsPane.setMinimumSize(new Dimension(550, 196));
                thresholdsPane.setLayout(null);

                //======== contentPanel ========
                {
                    contentPanel.setBorder(new EtchedBorder());
                    contentPanel.setLayout(null);

                    //---- label2 ----
                    label2.setText("Amplification:");
                    label2.setHorizontalAlignment(SwingConstants.RIGHT);
                    label2.setLabelFor(ampInput);
                    label2.setToolTipText("Amplification score, on a log-normalized scale");
                    label2.setPreferredSize(new Dimension(90, 18));
                    contentPanel.add(label2);
                    label2.setBounds(140, 96, label2.getPreferredSize().width, 18);

                    //---- label3 ----
                    label3.setText("Deletion:");
                    label3.setHorizontalAlignment(SwingConstants.RIGHT);
                    label3.setLabelFor(delInput);
                    label3.setToolTipText("Deletion score, on a log-normalized scale");
                    label3.setPreferredSize(new Dimension(60, 16));
                    contentPanel.add(label3);
                    label3.setBounds(360, 96, label3.getPreferredSize().width, 18);

                    //---- delInput ----
                    delInput.setText("0.9");
                    delInput.setMinimumSize(new Dimension(34, 28));
                    delInput.setPreferredSize(new Dimension(45, 28));
                    delInput.setMaximumSize(new Dimension(50, 2147483647));
                    contentPanel.add(delInput);
                    delInput.setBounds(new Rectangle(new Point(430, 91), delInput.getPreferredSize()));

                    //---- label4 ----
                    label4.setText("Up:");
                    label4.setHorizontalAlignment(SwingConstants.RIGHT);
                    label4.setLabelFor(expUpInput);
                    label4.setToolTipText("Expression score, log-normalized scale");
                    label4.setPreferredSize(new Dimension(100, 18));
                    contentPanel.add(label4);
                    label4.setBounds(130, 168, label4.getPreferredSize().width, 18);

                    //---- expUpInput ----
                    expUpInput.setText("1.0");
                    expUpInput.setMinimumSize(new Dimension(34, 28));
                    expUpInput.setPreferredSize(new Dimension(45, 28));
                    contentPanel.add(expUpInput);
                    expUpInput.setBounds(new Rectangle(new Point(240, 162), expUpInput.getPreferredSize()));

                    //---- label7 ----
                    label7.setText("Down:");
                    label7.setHorizontalAlignment(SwingConstants.RIGHT);
                    label7.setLabelFor(expDownInput);
                    label7.setToolTipText("Expression score, log-normalized scale");
                    label7.setPreferredSize(new Dimension(120, 16));
                    contentPanel.add(label7);
                    label7.setBounds(300, 168, label7.getPreferredSize().width, 18);

                    //---- expDownInput ----
                    expDownInput.setText("1.0");
                    expDownInput.setPreferredSize(new Dimension(45, 28));
                    expDownInput.setMinimumSize(new Dimension(34, 28));
                    expDownInput.setMaximumSize(new Dimension(50, 2147483647));
                    contentPanel.add(expDownInput);
                    expDownInput.setBounds(new Rectangle(new Point(430, 162), expDownInput.getPreferredSize()));

                    //---- ampInput ----
                    ampInput.setText("0.9");
                    ampInput.setMinimumSize(new Dimension(34, 28));
                    ampInput.setPreferredSize(new Dimension(45, 28));
                    contentPanel.add(ampInput);
                    ampInput.setBounds(new Rectangle(new Point(240, 91), ampInput.getPreferredSize()));

                    //---- label1 ----
                    label1.setText("Mutation:");
                    label1.setHorizontalAlignment(SwingConstants.RIGHT);
                    label1.setLabelFor(mutInput);
                    label1.setToolTipText("Minimum number of mutations found");
                    label1.setPreferredSize(new Dimension(66, 18));
                    contentPanel.add(label1);
                    label1.setBounds(50, 26, label1.getPreferredSize().width, 18);

                    //---- mutInput ----
                    mutInput.setText("1");
                    mutInput.setAutoscrolls(false);
                    mutInput.setMinimumSize(new Dimension(34, 28));
                    mutInput.setPreferredSize(new Dimension(45, 28));
                    contentPanel.add(mutInput);
                    mutInput.setBounds(new Rectangle(new Point(240, 21), mutInput.getPreferredSize()));

                    //---- label6 ----
                    label6.setText("Copy Number:");
                    contentPanel.add(label6);
                    label6.setBounds(30, 96, label6.getPreferredSize().width, 18);

                    //---- label8 ----
                    label8.setText("Expression:");
                    label8.setHorizontalAlignment(SwingConstants.RIGHT);
                    contentPanel.add(label8);
                    label8.setBounds(30, 168, 86, 18);
                    contentPanel.add(separator1);
                    separator1.setBounds(0, 135, 500, 10);
                    contentPanel.add(separator3);
                    separator3.setBounds(0, 65, 500, 10);

                    //---- separator2 ----
                    separator2.setPreferredSize(new Dimension(10, 210));
                    separator2.setOrientation(SwingConstants.VERTICAL);
                    contentPanel.add(separator2);
                    separator2.setBounds(new Rectangle(new Point(120, 0), separator2.getPreferredSize()));

                    { // compute preferred size
                        Dimension preferredSize = new Dimension();
                        for (int i = 0; i < contentPanel.getComponentCount(); i++) {
                            Rectangle bounds = contentPanel.getComponent(i).getBounds();
                            preferredSize.width = Math.max(bounds.x + bounds.width, preferredSize.width);
                            preferredSize.height = Math.max(bounds.y + bounds.height, preferredSize.height);
                        }
                        Insets insets = contentPanel.getInsets();
                        preferredSize.width += insets.right;
                        preferredSize.height += insets.bottom;
                        contentPanel.setMinimumSize(preferredSize);
                        contentPanel.setPreferredSize(preferredSize);
                    }
                }
                thresholdsPane.add(contentPanel);
                contentPanel.setBounds(12, 80, 500, 210);

                //======== panel2 ========
                {
                    panel2.setLayout(null);

                    { // compute preferred size
                        Dimension preferredSize = new Dimension();
                        for (int i = 0; i < panel2.getComponentCount(); i++) {
                            Rectangle bounds = panel2.getComponent(i).getBounds();
                            preferredSize.width = Math.max(bounds.x + bounds.width, preferredSize.width);
                            preferredSize.height = Math.max(bounds.y + bounds.height, preferredSize.height);
                        }
                        Insets insets = panel2.getInsets();
                        preferredSize.width += insets.right;
                        preferredSize.height += insets.bottom;
                        panel2.setMinimumSize(preferredSize);
                        panel2.setPreferredSize(preferredSize);
                    }
                }
                thresholdsPane.add(panel2);
                panel2.setBounds(new Rectangle(new Point(55, 25), panel2.getPreferredSize()));

                //---- textArea1 ----
                textArea1.setText("Samples are considered to have a given \"event\" if the value is above the thresholds below.");
                textArea1.setEditable(false);
                textArea1.setLineWrap(true);
                textArea1.setBackground(UIManager.getColor("Button.background"));
                thresholdsPane.add(textArea1);
                textArea1.setBounds(15, 10, 430, 40);

                //---- resetToDefaultsButton ----
                resetToDefaultsButton.setText("Reset to Defaults");
                resetToDefaultsButton.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        resetToDefaultsButtonActionPerformed(e);
                    }
                });
                thresholdsPane.add(resetToDefaultsButton);
                resetToDefaultsButton.setBounds(new Rectangle(new Point(10, 50), resetToDefaultsButton.getPreferredSize()));

                { // compute preferred size
                    Dimension preferredSize = new Dimension();
                    for (int i = 0; i < thresholdsPane.getComponentCount(); i++) {
                        Rectangle bounds = thresholdsPane.getComponent(i).getBounds();
                        preferredSize.width = Math.max(bounds.x + bounds.width, preferredSize.width);
                        preferredSize.height = Math.max(bounds.y + bounds.height, preferredSize.height);
                    }
                    Insets insets = thresholdsPane.getInsets();
                    preferredSize.width += insets.right;
                    preferredSize.height += insets.bottom;
                    thresholdsPane.setMinimumSize(preferredSize);
                    thresholdsPane.setPreferredSize(preferredSize);
                }
            }
            tabbedPane.addTab("Thresholds", thresholdsPane);

        }
        contentPane2.add(tabbedPane, BorderLayout.NORTH);
        pack();
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JTabbedPane tabbedPane;
    private JPanel seedGenesPane;
    private JPanel vSpacer1;
    private JLabel label5;
    private JScrollPane scrollPane2;
    private JTextArea seedGenesText;
    private JTextPane textPane1;
    private JPanel seedButtonBar;
    private JButton retrieveNetworkButton;
    private JButton resetSeedGeneTextAreaButton;
    private JButton cancel2;
    private JPanel filterPane;
    private JPanel panel1;
    private JButton addRow;
    private JPanel contentPane;
    private JScrollPane scrollPane1;
    private JTable geneTable;
    private JPanel buttonBar;
    private JLabel totNumGenes;
    private JCheckBox showSeedOnly;
    private JButton refFilter;
    private JCheckBox keepIsolated;
    private JButton okButton;
    private JButton cancelButton;
    private JButton saveButton;
    private JButton helpButton;
    private JPanel thresholdsPane;
    private JPanel contentPanel;
    private JLabel label2;
    private JLabel label3;
    private JTextField delInput;
    private JLabel label4;
    private JTextField expUpInput;
    private JLabel label7;
    private JTextField expDownInput;
    private JTextField ampInput;
    private JLabel label1;
    private JTextField mutInput;
    private JLabel label6;
    private JLabel label8;
    private JSeparator separator1;
    private JSeparator separator3;
    private JSeparator separator2;
    private JPanel panel2;
    private JTextArea textArea1;
    private JButton resetToDefaultsButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables


    /**
     * Export the current table to a tab-delimited file.
     * String exported should be the same as what user sees
     *
     * @param outFile
     * @throws IOException
     */
    private void saveTable(File outFile) throws FileNotFoundException {
        PrintWriter writer = new PrintWriter(outFile);
        TableModel model = geneTable.getModel();
        String delimiter = "\t";

        //Write header
        String header = model.getColumnName(0);
        for (int col = 1; col < model.getColumnCount(); col++) {
            header += delimiter + model.getColumnName(col);
        }

        writer.println(header);

        for (int row = 0; row < model.getRowCount(); row++) {
            String rowStr = "" + model.getValueAt(row, 0);
            for (int col = 1; col < model.getColumnCount(); col++) {
                rowStr += delimiter + model.getValueAt(row, col);
            }

            writer.println(rowStr);
        }

        writer.flush();
        writer.close();
    }

    private class GraphListModel extends AbstractTableModel {

        private List<Node> geneVertices = null;

        private List<Node> getGeneVertices() {
            if (geneVertices == null) {
                Collection<Node> nodes = network.geneVertexes();
                geneVertices = Arrays.asList(nodes.toArray(new Node[0]));
            }
            return geneVertices;
        }

        public void markDirty() {
            this.geneVertices = null;
            this.fireTableStructureChanged();
        }

        @Override
        public int getRowCount() {
            return getGeneVertices().size();
        }

        @Override
        public Class<?> getColumnClass(int columnIndex) {
            switch (columnIndex) {
                case 0:
                    return String.class;
                case 1:
                    return Integer.class;
                default:
                    return Double.class;
            }
        }

        @Override
        public int getColumnCount() {
            return columnNames.size();
        }

        public String getColumnName(int col) {
            return columnNames.get(col);
        }

        @Override
        public Object getValueAt(int rowIndex, int columnIndex) {

            Node n = getGeneVertices().get(rowIndex);
            String nm = GeneNetwork.getNodeKeyData(n, "label");
            switch (columnIndex) {
                case 0:
                    return nm;
                case 1:
                    return network.edgesOf(n).size();
                default:
                    String key = columnNumToKeyMap.get(columnIndex);
                    if (key == null) {
                        return null;
                    }
                    String val = GeneNetwork.getNodeKeyData(n, key);
                    if ("nan".equalsIgnoreCase(val) || val == null) {
                        return null;
                    }

                    //Change from fraction to percent
                    double dPerc = Double.parseDouble(val) * 100;
                    return dPerc;
            }
        }

        @Override
        public boolean isCellEditable(int row, int col) {
            return false;
        }
    }

    private class RefreshListener implements FocusListener, ActionListener {

        @Override
        public void focusGained(FocusEvent e) {
            //pass
        }

        @Override
        public void focusLost(FocusEvent e) {
            refreshFilters();
        }

        @Override
        public void actionPerformed(ActionEvent e) {
            refreshFilters();
        }
    }
}
