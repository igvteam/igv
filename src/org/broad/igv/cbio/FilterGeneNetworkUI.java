/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

/*
 * Created by JFormDesigner on Tue Mar 06 14:09:15 EST 2012
 */

package org.broad.igv.cbio;

import org.apache.commons.collections.Predicate;
import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.PreferenceManager;
import org.broad.igv.lists.GeneList;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.ui.util.*;
import org.broad.igv.util.BrowserLauncher;
import org.broad.igv.util.HttpUtils;
import org.w3c.dom.Node;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.border.EtchedBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
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

    private static final int MAX_GENES_JUSTSHOW = 20;

    private GeneList geneList;
    private List<AttributeFilter> filterRows = new ArrayList<AttributeFilter>(1);
    GeneNetwork network = null;

    private GraphListModel listModel;
    private Map<String, JTextField> thresholdsMap = new HashMap<String, JTextField>(5);


    private static List<String> columnNames;
    private static Map<Integer, String> columnNumToKeyMap;

    static {
        String[] firstLabels = {"Gene label", "Interactions"};
        columnNumToKeyMap = new HashMap<Integer, String>(GeneNetwork.attributeMap.size());
        columnNames = new ArrayList<String>(firstLabels.length + GeneNetwork.attributeMap.size());
        for (String label : firstLabels) {
            columnNames.add(label);
        }
        int ind = columnNames.size();
        for (String label : GeneNetwork.attributeMap.keySet()) {
            columnNumToKeyMap.put(ind, label);
            ind++;

            label = label.replace('_', ' ');
            label = label.replace("PERCENT", "%");
            columnNames.add(label);

        }

    }

    public FilterGeneNetworkUI(Frame owner, GeneList geneList) {
        super(owner);
        this.geneList = geneList;
    }

    @Override
    public void setVisible(boolean visible) {
        if (visible && network == null) {
            loadcBioData(this.geneList.getLoci());

        } else {
            //Just retrieved network; if
            //there are only a small number of genes
            //we skip to it
            if (network.vertexSet().size() < MAX_GENES_JUSTSHOW) {
                showNetwork();
            } else {
                super.setVisible(visible);
            }
        }
    }

    private void loadcBioData(final List<String> geneLoci) {

        final IndefiniteProgressMonitor indefMonitor = new IndefiniteProgressMonitor(60);
        final ProgressBar progressBar = ProgressBar.showProgressDialog((Frame) getOwner(), "Loading cBio data...", indefMonitor, true);
        progressBar.setIndeterminate(true);
        indefMonitor.start();

        final Runnable showUI = new Runnable() {
            @Override
            public void run() {
                initComponents();
                initComponentData();
                setVisible(true);
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
                        MessageUtils.showMessage("No results found for " + HttpUtils.buildURLString(geneLoci, ", "));
                    } else {
                        network.annotateAll(IGV.getInstance().getAllTracks(false));
                        UIUtilities.invokeOnEventThread(showUI);
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                    log.error(e.getMessage());
                    MessageUtils.showMessage("Error loading data: " + e.getMessage());
                } finally {
                    WaitCursorManager.removeWaitCursor(token);

                    if (progressBar != null) {
                        progressBar.close();
                        indefMonitor.stop();
                    }
                }
            }
        };

        // If we're on the dispatch thread spawn a worker, otherwise just execute.
        if (SwingUtilities.isEventDispatchThread()) {
            SwingWorker worker = new SwingWorker() {
                @Override
                protected Object doInBackground() throws Exception {
                    runnable.run();
                    return null;
                }
            };

            worker.execute();
        } else {
            runnable.run();
        }
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
        add();

        initThresholdsMap();

        loadThresholds();

        listModel = new GraphListModel();
        geneTable.setModel(listModel);
        applySoftFilters();
    }

    private void remove(AttributeFilter row) {
        contentPane.remove(row.getPanel());
        filterRows.remove(row);
        int numRows = filterRows.size();
        filterRows.get(numRows - 1).setIsLast(true);
        filterRows.get(0).setShowDel(numRows >= 2);
        validateTree();
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

        validateTree();

    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        setVisible(false);
    }

    private void applySoftFilters() {
        network.clearAllFilters();

        //TODO This is only AND, should also include OR
        for (AttributeFilter filter : this.filterRows) {
            String filt_el = (String) filter.getAttrName().getSelectedItem();
            if (GeneNetwork.attributeMap.containsKey(filt_el) || GeneNetwork.PERCENT_ALTERED.equals(filt_el)) {
                float min = Float.parseFloat(filter.minVal.getText());
                float max = Float.parseFloat(filter.maxVal.getText());
                network.filterNodesRange(filt_el, min / 100, max / 100);
            }
        }
        if (!keepIsolated.isSelected()) {
            network.pruneGraph();
        }

        totNumGenes.setText("Total Genes: " + network.vertexSetFiltered().size());

        this.listModel.markDirty();
    }

    /**
     * TODO This should run on a separate thread
     */
    private void showNetwork() {

        try {
            String url = network.outputForcBioView();
            url = "file://" + url;
            BrowserLauncher.openURL(url);
        } catch (IOException err) {
            MessageUtils.showMessage("Error opening network for viewing. " + err.getMessage());
        }
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
        //If any rows are selected, we only keep those
        //We lookup by name, because table could get sorted
        int[] keepRows = getModelIndices();
        if (keepRows.length > 0) {
            final Set<Node> keepNodes = new HashSet<Node>(keepRows.length);
            GraphListModel model = (GraphListModel) geneTable.getModel();
            List<Node> vertices = model.getVertices();
            for (Integer loc : keepRows) {
                keepNodes.add(vertices.get(loc));
            }

            Predicate<Node> selectedPredicated = new Predicate<Node>() {

                @Override
                public boolean evaluate(Node object) {
                    return keepNodes.contains(object);
                }
            };
            network.filterNodes(selectedPredicated);
        }

        setVisible(false);
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
        this.validateTree();
    }

    private void refFilterActionPerformed(ActionEvent e) {
        refreshFilters();
    }

    private boolean saveThresholds() {
        try {
            for (Map.Entry<String, JTextField> entry : thresholdsMap.entrySet()) {
                float value = Float.parseFloat(entry.getValue().getText());
                PreferenceManager.getInstance().put(entry.getKey(), "" + value);
            }
        } catch (NumberFormatException e) {
            MessageUtils.showMessage("Inputs must be numeric. " + e.getMessage());
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
        int thresholdTabNum = tabbedPane.indexOfTab("Thresholds");
        if (thresholdTabNum < 0) {
            //Component not built yet
            return;
        }
        Component thresholdTab = tabbedPane.getComponentAt(thresholdTabNum);
        if (!tabbedPane.getSelectedComponent().equals(thresholdTab) && !saveThresholds()) {
            tabbedPane.setSelectedComponent(thresholdTab);
        }
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

    private void exportButtonActionPerformed(ActionEvent e) {
        // TODO add your code here
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        tabbedPane = new JTabbedPane();
        dialogPane = new JPanel();
        panel1 = new JPanel();
        addRow = new JButton();
        contentPane = new JPanel();
        scrollPane1 = new JScrollPane();
        geneTable = new JTable();
        panel3 = new JPanel();
        buttonBar = new JPanel();
        totNumGenes = new JLabel();
        keepIsolated = new JCheckBox();
        okButton = new JButton();
        refFilter = new JButton();
        cancelButton = new JButton();
        helpButton = new JButton();
        saveButton = new JButton();
        thresholds = new JPanel();
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
        panel2 = new JPanel();
        textArea1 = new JTextArea();

        //======== this ========
        setMinimumSize(new Dimension(600, 22));
        setModalityType(Dialog.ModalityType.DOCUMENT_MODAL);
        Container contentPane2 = getContentPane();
        contentPane2.setLayout(new BorderLayout());

        //======== tabbedPane ========
        {
            tabbedPane.setPreferredSize(new Dimension(550, 346));
            tabbedPane.setMinimumSize(new Dimension(550, 346));
            tabbedPane.addChangeListener(new ChangeListener() {
                @Override
                public void stateChanged(ChangeEvent e) {
                    tabbedPaneStateChanged(e);
                }
            });

            //======== dialogPane ========
            {
                dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
                dialogPane.setMinimumSize(new Dimension(443, 300));
                dialogPane.setPreferredSize(new Dimension(443, 300));
                dialogPane.setLayout(new GridBagLayout());
                ((GridBagLayout) dialogPane.getLayout()).columnWidths = new int[]{0, 0};
                ((GridBagLayout) dialogPane.getLayout()).rowHeights = new int[]{0, 0, 0, 0, 0, 0};
                ((GridBagLayout) dialogPane.getLayout()).columnWeights = new double[]{1.0, 1.0E-4};
                ((GridBagLayout) dialogPane.getLayout()).rowWeights = new double[]{0.0, 0.0, 1.0, 0.0, 0.0, 1.0E-4};

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
                dialogPane.add(panel1, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 0, 0), 0, 0));

                //======== contentPane ========
                {
                    contentPane.setLayout(new BoxLayout(contentPane, BoxLayout.Y_AXIS));
                }
                dialogPane.add(contentPane, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 0, 0), 0, 0));

                //======== scrollPane1 ========
                {

                    //---- geneTable ----
                    geneTable.setAutoCreateRowSorter(true);
                    scrollPane1.setViewportView(geneTable);
                }
                dialogPane.add(scrollPane1, new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 0, 0), 0, 0));

                //======== panel3 ========
                {
                    panel3.setLayout(null);

                    { // compute preferred size
                        Dimension preferredSize = new Dimension();
                        for (int i = 0; i < panel3.getComponentCount(); i++) {
                            Rectangle bounds = panel3.getComponent(i).getBounds();
                            preferredSize.width = Math.max(bounds.x + bounds.width, preferredSize.width);
                            preferredSize.height = Math.max(bounds.y + bounds.height, preferredSize.height);
                        }
                        Insets insets = panel3.getInsets();
                        preferredSize.width += insets.right;
                        preferredSize.height += insets.bottom;
                        panel3.setMinimumSize(preferredSize);
                        panel3.setPreferredSize(preferredSize);
                    }
                }
                dialogPane.add(panel3, new GridBagConstraints(0, 4, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 0, 0), 0, 0));

                //======== buttonBar ========
                {
                    buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                    buttonBar.setLayout(new GridBagLayout());
                    ((GridBagLayout) buttonBar.getLayout()).columnWidths = new int[]{0, 85, 85, 80};
                    ((GridBagLayout) buttonBar.getLayout()).columnWeights = new double[]{1.0, 0.0, 0.0, 0.0};

                    //---- totNumGenes ----
                    totNumGenes.setText("Total Genes: #");
                    buttonBar.add(totNumGenes, new GridBagConstraints(0, 0, 2, 1, 0.0, 0.0,
                            GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                            new Insets(0, 0, 5, 5), 0, 0));

                    //---- keepIsolated ----
                    keepIsolated.setText("Keep Isolated Genes");
                    keepIsolated.setVisible(false);
                    buttonBar.add(keepIsolated, new GridBagConstraints(0, 3, 2, 1, 0.0, 0.0,
                            GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                            new Insets(0, 0, 5, 5), 0, 0));

                    //---- okButton ----
                    okButton.setText("View Network");
                    okButton.setToolTipText("Display the network in a web browser");
                    okButton.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            okButtonActionPerformed(e);
                        }
                    });
                    buttonBar.add(okButton, new GridBagConstraints(0, 4, 1, 1, 0.0, 0.0,
                            GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                            new Insets(0, 0, 0, 5), 0, 0));

                    //---- refFilter ----
                    refFilter.setText("Refresh Filter");
                    refFilter.setVisible(false);
                    refFilter.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            refFilterActionPerformed(e);
                        }
                    });
                    buttonBar.add(refFilter, new GridBagConstraints(1, 2, 1, 1, 0.0, 0.0,
                            GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                            new Insets(0, 0, 5, 5), 0, 0));

                    //---- cancelButton ----
                    cancelButton.setText("Cancel");
                    cancelButton.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            cancelButtonActionPerformed(e);
                        }
                    });
                    buttonBar.add(cancelButton, new GridBagConstraints(3, 4, 1, 1, 0.0, 0.0,
                            GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                            new Insets(0, 0, 0, 0), 0, 0));

                    //---- helpButton ----
                    helpButton.setText("Help");
                    helpButton.setVisible(false);
                    buttonBar.add(helpButton, new GridBagConstraints(3, 2, 1, 1, 0.0, 0.0,
                            GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                            new Insets(0, 0, 5, 0), 0, 0));

                    //---- saveButton ----
                    saveButton.setText("Save Table");
                    saveButton.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            saveButtonActionPerformed(e);
                        }
                    });
                    buttonBar.add(saveButton, new GridBagConstraints(1, 4, 1, 1, 0.0, 0.0,
                            GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                            new Insets(0, 0, 0, 5), 0, 0));
                }
                dialogPane.add(buttonBar, new GridBagConstraints(0, 3, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 0, 0), 0, 0));
            }
            tabbedPane.addTab("Filter", dialogPane);


            //======== thresholds ========
            {
                thresholds.setPreferredSize(new Dimension(550, 196));
                thresholds.setMinimumSize(new Dimension(550, 196));
                thresholds.setLayout(null);

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
                    delInput.setText("0.7");
                    delInput.setMinimumSize(new Dimension(34, 28));
                    delInput.setPreferredSize(new Dimension(45, 28));
                    delInput.setMaximumSize(new Dimension(50, 2147483647));
                    contentPanel.add(delInput);
                    delInput.setBounds(new Rectangle(new Point(240, 162), delInput.getPreferredSize()));

                    //---- label4 ----
                    label4.setText("Up:");
                    label4.setHorizontalAlignment(SwingConstants.RIGHT);
                    label4.setLabelFor(expUpInput);
                    label4.setToolTipText("Expression score, log-normalized scale");
                    label4.setPreferredSize(new Dimension(100, 18));
                    contentPanel.add(label4);
                    label4.setBounds(130, 168, label4.getPreferredSize().width, 18);

                    //---- expUpInput ----
                    expUpInput.setText("0.1");
                    expUpInput.setMinimumSize(new Dimension(34, 28));
                    expUpInput.setPreferredSize(new Dimension(45, 28));
                    contentPanel.add(expUpInput);
                    expUpInput.setBounds(new Rectangle(new Point(430, 91), expUpInput.getPreferredSize()));

                    //---- label7 ----
                    label7.setText("Down:");
                    label7.setHorizontalAlignment(SwingConstants.RIGHT);
                    label7.setLabelFor(expDownInput);
                    label7.setToolTipText("Expression score, log-normalized scale");
                    label7.setPreferredSize(new Dimension(120, 16));
                    contentPanel.add(label7);
                    label7.setBounds(300, 168, label7.getPreferredSize().width, 18);

                    //---- expDownInput ----
                    expDownInput.setText("0.1");
                    expDownInput.setPreferredSize(new Dimension(45, 28));
                    expDownInput.setMinimumSize(new Dimension(34, 28));
                    expDownInput.setMaximumSize(new Dimension(50, 2147483647));
                    contentPanel.add(expDownInput);
                    expDownInput.setBounds(new Rectangle(new Point(430, 162), expDownInput.getPreferredSize()));

                    //---- ampInput ----
                    ampInput.setText("0.7");
                    ampInput.setMinimumSize(new Dimension(34, 28));
                    ampInput.setPreferredSize(new Dimension(45, 28));
                    contentPanel.add(ampInput);
                    ampInput.setBounds(new Rectangle(new Point(240, 91), ampInput.getPreferredSize()));

                    //---- label1 ----
                    label1.setText("Mutation:");
                    label1.setHorizontalAlignment(SwingConstants.RIGHT);
                    label1.setLabelFor(mutInput);
                    label1.setToolTipText("Minimum number of mutations found");
                    label1.setPreferredSize(new Dimension(65, 18));
                    contentPanel.add(label1);
                    label1.setBounds(165, 26, label1.getPreferredSize().width, 18);

                    //---- mutInput ----
                    mutInput.setText("1");
                    mutInput.setAutoscrolls(false);
                    mutInput.setMinimumSize(new Dimension(34, 28));
                    mutInput.setPreferredSize(new Dimension(45, 28));
                    contentPanel.add(mutInput);
                    mutInput.setBounds(new Rectangle(new Point(240, 21), mutInput.getPreferredSize()));

                    //---- label6 ----
                    label6.setText("Copy Number");
                    contentPanel.add(label6);
                    label6.setBounds(30, 96, label6.getPreferredSize().width, 18);

                    //---- label8 ----
                    label8.setText("Expression");
                    label8.setHorizontalAlignment(SwingConstants.RIGHT);
                    contentPanel.add(label8);
                    label8.setBounds(30, 168, 86, 18);
                    contentPanel.add(separator1);
                    separator1.setBounds(0, 135, 500, 10);
                    contentPanel.add(separator3);
                    separator3.setBounds(0, 65, 500, 10);

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
                thresholds.add(contentPanel);
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
                thresholds.add(panel2);
                panel2.setBounds(new Rectangle(new Point(55, 25), panel2.getPreferredSize()));

                //---- textArea1 ----
                textArea1.setText("Samples are considered to have a given \"event\" if the value is above the thresholds below.");
                textArea1.setEditable(false);
                textArea1.setLineWrap(true);
                textArea1.setBackground(UIManager.getColor("Button.background"));
                thresholds.add(textArea1);
                textArea1.setBounds(15, 10, 430, 40);

                { // compute preferred size
                    Dimension preferredSize = new Dimension();
                    for (int i = 0; i < thresholds.getComponentCount(); i++) {
                        Rectangle bounds = thresholds.getComponent(i).getBounds();
                        preferredSize.width = Math.max(bounds.x + bounds.width, preferredSize.width);
                        preferredSize.height = Math.max(bounds.y + bounds.height, preferredSize.height);
                    }
                    Insets insets = thresholds.getInsets();
                    preferredSize.width += insets.right;
                    preferredSize.height += insets.bottom;
                    thresholds.setMinimumSize(preferredSize);
                    thresholds.setPreferredSize(preferredSize);
                }
            }
            tabbedPane.addTab("Thresholds", thresholds);

        }
        contentPane2.add(tabbedPane, BorderLayout.NORTH);
        pack();
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JTabbedPane tabbedPane;
    private JPanel dialogPane;
    private JPanel panel1;
    private JButton addRow;
    private JPanel contentPane;
    private JScrollPane scrollPane1;
    private JTable geneTable;
    private JPanel panel3;
    private JPanel buttonBar;
    private JLabel totNumGenes;
    private JCheckBox keepIsolated;
    private JButton okButton;
    private JButton refFilter;
    private JButton cancelButton;
    private JButton helpButton;
    private JButton saveButton;
    private JPanel thresholds;
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
    private JPanel panel2;
    private JTextArea textArea1;
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

        private List<Node> vertices = null;

        private List<Node> getVertices() {
            if (vertices == null) {
                Set<Node> nodes = network.vertexSetFiltered();
                vertices = Arrays.asList(nodes.toArray(new Node[0]));
            }
            return vertices;

            //Collections.sort(vertexNames);
        }

        public void markDirty() {
            this.vertices = null;
            this.fireTableStructureChanged();
        }

        @Override
        public int getRowCount() {
            return getVertices().size();
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

            Node n = getVertices().get(rowIndex);
            String nm = GeneNetwork.getNodeKeyData(n, "label");
            switch (columnIndex) {
                case 0:
                    return nm;
                case 1:
                    return network.edgesOfFiltered(n).size();
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

                    if (dPerc == 0.0d) return "0.0";
                    //If above 1, just show integer. If small, show in exponential format
                    String fmt = "%2.1f";
                    if (dPerc < 0.1d) {
                        fmt = "%2.1e";
                    }

                    String sVal = String.format(fmt, dPerc);
                    return sVal;

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
