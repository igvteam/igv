/*
 * Created by JFormDesigner on Tue Mar 06 14:09:15 EST 2012
 */

package org.broad.igv.cbio;

import org.broad.igv.PreferenceManager;
import org.broad.igv.lists.GeneList;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.BrowserLauncher;
import org.w3c.dom.Node;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

/**
 * Dialog for letting the user filter a GeneNetwork from
 * cBio.
 *
 * @author Jacob Silterra
 */
public class FilterGeneNetworkUI extends JDialog {

    private GeneList geneList;
    private List<AttributeFilter> filterRows = new ArrayList<AttributeFilter>(1);
    GeneNetwork network;

    private GraphListModel listModel;

    public FilterGeneNetworkUI(Frame owner, GeneList geneList) {
        super(owner);
        this.geneList = geneList;
        initComponents();
        initComponentData(geneList.getLoci());
    }

    /**
     * Load relevant data from IGV to initialize
     * displayed components.
     */
    private void initComponentData(List<String> geneLoci) {
        add();
        loadThresholds();

        network = GeneNetwork.getFromCBIO(geneLoci);
        network.annotateAll(IGV.getInstance().getAllTracks(false));
        listModel = new GraphListModel();
        geneTable.setModel(listModel);
        applySoftFilters();
    }

    private void add() {
        final AttributeFilter row = new AttributeFilter();

        row.getDelRow().addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                contentPane.remove(row.getComponent());
                filterRows.remove(row);
                validateTree();
            }
        });
        contentPane.add(row.getComponent());
        validateTree();

        this.filterRows.add(row);
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        setVisible(false);
    }

    private void applySoftFilters() {
        network.clearAllFilters();

        //TODO This is only AND, should also include OR
        for (AttributeFilter filter : this.filterRows) {
            String filt_el = (String) filter.attrName.getSelectedItem();
            if (GeneNetwork.attributeMap.containsKey(filt_el) || GeneNetwork.PERCENT_ALTERED.equals(filt_el)) {
                float min = Float.parseFloat(filter.minVal.getText());
                float max = Float.parseFloat(filter.maxVal.getText());
                network.filterNodesRange(filt_el, min / 100, max / 100);
            }
        }
        if (!keepIsolated.isSelected()) {
            network.pruneGraph();
        }

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

    private void okButtonActionPerformed(ActionEvent e) {
        setVisible(false);
        showNetwork();
    }

    private void addRowActionPerformed(ActionEvent e) {
        add();
    }

    private void refFilterActionPerformed(ActionEvent e) {
        this.applySoftFilters();
        listModel.markDirty();
        this.validateTree();

    }

    private boolean saveThresholds() {
        try {
            float mut = Float.parseFloat(mutInput.getText());
            float amp = Float.parseFloat(ampInput.getText());
            float exp = Float.parseFloat(expInput.getText());
            PreferenceManager.getInstance().put(PreferenceManager.CBIO_MUTATION_THRESHOLD, "" + mut);
            PreferenceManager.getInstance().put(PreferenceManager.CBIO_AMPLIFICATION_THRESHOLD, "" + amp);
            PreferenceManager.getInstance().put(PreferenceManager.CBIO_EXPRESSION_THRESHOLD, "" + exp);
        } catch (NumberFormatException e) {
            MessageUtils.showMessage("Inputs must be numeric. " + e.getMessage());
            return false;
        }
        return true;
    }

    private void loadThresholds() {
        String mut = PreferenceManager.getInstance().get(PreferenceManager.CBIO_MUTATION_THRESHOLD, "" + 0.1);
        String amp = PreferenceManager.getInstance().get(PreferenceManager.CBIO_AMPLIFICATION_THRESHOLD, "" + 0.1);
        String exp = PreferenceManager.getInstance().get(PreferenceManager.CBIO_EXPRESSION_THRESHOLD, "" + 0.1);
        mutInput.setText(mut);
        ampInput.setText(amp);
        expInput.setText(exp);
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

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        tabbedPane = new JTabbedPane();
        dialogPane = new JPanel();
        contentPane = new JPanel();
        buttonBar = new JPanel();
        addRow = new JButton();
        keepIsolated = new JCheckBox();
        okButton = new JButton();
        refFilter = new JButton();
        cancelButton = new JButton();
        helpButton = new JButton();
        scrollPane1 = new JScrollPane();
        geneTable = new JTable();
        thresholds = new JPanel();
        contentPanel = new JPanel();
        label5 = new JLabel();
        label6 = new JLabel();
        label1 = new JLabel();
        mutInput = new JTextField();
        label2 = new JLabel();
        ampInput = new JTextField();
        label4 = new JLabel();
        expInput = new JTextField();

        //======== this ========
        setMinimumSize(new Dimension(550, 22));
        Container contentPane2 = getContentPane();
        contentPane2.setLayout(new BorderLayout());

        //======== tabbedPane ========
        {
            tabbedPane.setPreferredSize(new Dimension(464, 346));
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
                ((GridBagLayout) dialogPane.getLayout()).rowHeights = new int[]{0, 0, 0, 0};
                ((GridBagLayout) dialogPane.getLayout()).columnWeights = new double[]{1.0, 1.0E-4};
                ((GridBagLayout) dialogPane.getLayout()).rowWeights = new double[]{0.0, 1.0, 0.0, 1.0E-4};

                //======== contentPane ========
                {
                    contentPane.setLayout(new BoxLayout(contentPane, BoxLayout.Y_AXIS));
                }
                dialogPane.add(contentPane, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 0, 0), 0, 0));

                //======== buttonBar ========
                {
                    buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                    buttonBar.setLayout(new GridBagLayout());
                    ((GridBagLayout) buttonBar.getLayout()).columnWidths = new int[]{0, 85, 85, 80};
                    ((GridBagLayout) buttonBar.getLayout()).columnWeights = new double[]{1.0, 0.0, 0.0, 0.0};

                    //---- addRow ----
                    addRow.setText("Add");
                    addRow.setMaximumSize(new Dimension(30, 28));
                    addRow.setMinimumSize(new Dimension(30, 28));
                    addRow.setPreferredSize(new Dimension(30, 28));
                    addRow.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            addRowActionPerformed(e);
                        }
                    });
                    buttonBar.add(addRow, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0,
                            GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                            new Insets(0, 0, 5, 5), 0, 0));

                    //---- keepIsolated ----
                    keepIsolated.setText("Keep Isolated Genes");
                    buttonBar.add(keepIsolated, new GridBagConstraints(0, 1, 2, 1, 0.0, 0.0,
                            GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                            new Insets(0, 0, 5, 5), 0, 0));

                    //---- okButton ----
                    okButton.setText("View");
                    okButton.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            okButtonActionPerformed(e);
                        }
                    });
                    buttonBar.add(okButton, new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0,
                            GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                            new Insets(0, 0, 0, 5), 0, 0));

                    //---- refFilter ----
                    refFilter.setText("Refresh Filter");
                    refFilter.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            refFilterActionPerformed(e);
                        }
                    });
                    buttonBar.add(refFilter, new GridBagConstraints(1, 2, 1, 1, 0.0, 0.0,
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
                    buttonBar.add(cancelButton, new GridBagConstraints(2, 2, 1, 1, 0.0, 0.0,
                            GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                            new Insets(0, 0, 0, 5), 0, 0));

                    //---- helpButton ----
                    helpButton.setText("Help");
                    helpButton.setVisible(false);
                    buttonBar.add(helpButton, new GridBagConstraints(3, 2, 1, 1, 0.0, 0.0,
                            GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                            new Insets(0, 0, 0, 0), 0, 0));
                }
                dialogPane.add(buttonBar, new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 0, 0), 0, 0));

                //======== scrollPane1 ========
                {
                    scrollPane1.setViewportView(geneTable);
                }
                dialogPane.add(scrollPane1, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 0, 0), 0, 0));
            }
            tabbedPane.addTab("Filter", dialogPane);


            //======== thresholds ========
            {
                thresholds.setBorder(new EmptyBorder(12, 12, 12, 12));
                thresholds.setLayout(new BorderLayout());

                //======== contentPanel ========
                {
                    contentPanel.setLayout(new GridLayout(4, 2, 20, 20));

                    //---- label5 ----
                    label5.setText("Minimum thresholds for a");
                    label5.setHorizontalAlignment(SwingConstants.RIGHT);
                    label5.setMaximumSize(new Dimension(361, 16));
                    contentPanel.add(label5);

                    //---- label6 ----
                    label6.setText("gene to be considered \"altered\".");
                    contentPanel.add(label6);

                    //---- label1 ----
                    label1.setText("Mutation:");
                    label1.setHorizontalAlignment(SwingConstants.RIGHT);
                    contentPanel.add(label1);

                    //---- mutInput ----
                    mutInput.setText("0.1");
                    contentPanel.add(mutInput);

                    //---- label2 ----
                    label2.setText("Amplification/Deletion:");
                    label2.setHorizontalAlignment(SwingConstants.RIGHT);
                    contentPanel.add(label2);

                    //---- ampInput ----
                    ampInput.setText("0.1");
                    contentPanel.add(ampInput);

                    //---- label4 ----
                    label4.setText("Expression Up/Down:");
                    label4.setHorizontalAlignment(SwingConstants.RIGHT);
                    contentPanel.add(label4);

                    //---- expInput ----
                    expInput.setText("0.1");
                    contentPanel.add(expInput);
                }
                thresholds.add(contentPanel, BorderLayout.CENTER);
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
    private JPanel contentPane;
    private JPanel buttonBar;
    private JButton addRow;
    private JCheckBox keepIsolated;
    private JButton okButton;
    private JButton refFilter;
    private JButton cancelButton;
    private JButton helpButton;
    private JScrollPane scrollPane1;
    private JTable geneTable;
    private JPanel thresholds;
    private JPanel contentPanel;
    private JLabel label5;
    private JLabel label6;
    private JLabel label1;
    private JTextField mutInput;
    private JLabel label2;
    private JTextField ampInput;
    private JLabel label4;
    private JTextField expInput;
    // JFormDesigner - End of variables declaration  //GEN-END:variables


    private class GraphListModel extends AbstractTableModel {

        private List<Node> vertices = null;
        private final String[] columnNames = {"Gene label", "# of interactions"};

        private void getVertices() {
            Set<Node> nodes = network.vertexSetFiltered();
            vertices = Arrays.asList(nodes.toArray(new Node[0]));

            //Collections.sort(vertexNames);
        }

        public void markDirty() {
            this.vertices = null;
            this.fireTableStructureChanged();
        }

        @Override
        public int getRowCount() {
            if (vertices == null) {
                getVertices();
            }
            return vertices.size();
        }

        @Override
        public int getColumnCount() {
            return 2;
        }

        public String getColumnName(int col) {
            return columnNames[col].toString();
        }

        @Override
        public Object getValueAt(int rowIndex, int columnIndex) {
            if (vertices == null) {
                getVertices();
            }

            Node n = vertices.get(rowIndex);
            String nm = GeneNetwork.getNodeKeyData(vertices.get(rowIndex), "label");
            switch (columnIndex) {
                case 0:
                    return nm;
                case 1:
                    return network.edgesOf(n).size();
                default:
                    return null;
            }
        }

        @Override
        public boolean isCellEditable(int row, int col) {
            return false;
        }
    }
}
