/*
 * Created by JFormDesigner on Tue Mar 06 14:09:15 EST 2012
 */

package org.broad.igv.cbio;

import org.broad.igv.lists.GeneList;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.BrowserLauncher;
import org.w3c.dom.Node;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
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


    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
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

        //======== this ========
        setMinimumSize(new Dimension(550, 22));
        Container contentPane2 = getContentPane();
        contentPane2.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setMinimumSize(new Dimension(443, 300));
            dialogPane.setPreferredSize(new Dimension(443, 300));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPane ========
            {
                contentPane.setLayout(new BoxLayout(contentPane, BoxLayout.Y_AXIS));
            }
            dialogPane.add(contentPane, BorderLayout.NORTH);

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
            dialogPane.add(buttonBar, BorderLayout.SOUTH);

            //======== scrollPane1 ========
            {
                scrollPane1.setViewportView(geneTable);
            }
            dialogPane.add(scrollPane1, BorderLayout.CENTER);
        }
        contentPane2.add(dialogPane, BorderLayout.CENTER);
        pack();
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
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
