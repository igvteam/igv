/*
 * Created by JFormDesigner on Tue Mar 06 14:09:15 EST 2012
 */

package org.broad.igv.cbio;

import org.broad.igv.lists.GeneList;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.BrowserLauncher;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Dialog for letting the user filter a GeneNetwork from
 * cBio.
 *
 * @author Jacob Silterra
 */
public class FilterGeneNetworkUI extends JDialog {

    private GeneList geneList;
    private List<AttributeFilter> filterRows = new ArrayList<AttributeFilter>(1);

    public FilterGeneNetworkUI(Frame owner, GeneList geneList) {
        super(owner);
        this.geneList = geneList;
        initComponents();
        initComponentData();
    }

    /**
     * Load relevant data from IGV to initialize
     * displayed components.
     */
    private void initComponentData() {
        add();
    }

    private void add() {
        final AttributeFilter row = new AttributeFilter();

        row.getDelRow().addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                contentPane.remove(row.getComponent());
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

    /**
     * TODO This should run on a separate thread
     *
     * @param geneLoci
     */
    private void showNetwork(List<String> geneLoci) {
        GeneNetwork network = GeneNetwork.getFromCBIO(geneLoci);

        network.annotateAll(IGV.getInstance().getAllTracks(false));

        //TODO This is only AND, should also include OR
        for (AttributeFilter filter : this.filterRows) {
            String filt_el = (String) filter.attrName.getSelectedItem();
            if (GeneNetwork.attribute_map.containsKey(filt_el) || GeneNetwork.PERCENT_ALTERED.equals(filt_el)) {
                float min = Float.parseFloat(filter.minVal.getText());
                float max = Float.parseFloat(filter.maxVal.getText());
                network.filterNodesRange(filt_el, min / 100, max / 100);
            }
        }
        if (!keepIsolated.isSelected()) {
            network.pruneGraph();
        }

        try {
            String url = network.outputForcBioView();
            url = "file://" + url;
            BrowserLauncher.openURL(url);
        } catch (IOException err) {
            MessageUtils.showMessage("Error opening network for viewing. " + err.getMessage());
        }
    }

    private void okButtonActionPerformed(ActionEvent e) {
        List<String> geneLoci = geneList.getLoci();
        setVisible(false);
        showNetwork(geneLoci);
    }

    private void addRowActionPerformed(ActionEvent e) {
        add();
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
        cancelButton = new JButton();
        helpButton = new JButton();

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
                buttonBar.add(keepIsolated, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 5, 5), 0, 0));

                //---- okButton ----
                okButton.setText("OK");
                okButton.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        okButtonActionPerformed(e);
                    }
                });
                buttonBar.add(okButton, new GridBagConstraints(1, 2, 1, 1, 0.0, 0.0,
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
    private JButton cancelButton;
    private JButton helpButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
