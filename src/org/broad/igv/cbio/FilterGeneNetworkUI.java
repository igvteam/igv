/*
 * Created by JFormDesigner on Tue Mar 06 14:09:15 EST 2012
 */

package org.broad.igv.cbio;

import com.jidesoft.swing.RangeSlider;
import org.broad.igv.lists.GeneList;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.BrowserLauncher;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.io.IOException;
import java.util.List;

/**
 * Dialog for letting the user filter a GeneNetwork from
 * cBio.
 *
 * @author Jacob Silterra
 */
public class FilterGeneNetworkUI extends JDialog {

    private GeneList geneList;

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
        for (String attrs : GeneNetwork.attribute_map.keySet()) {
            comboBox1.addItem(attrs);
        }
        comboBox1.addItem("None");
    }


    private void comboBox1ActionPerformed(ActionEvent e) {
        // TODO add your code here
    }

    private void comboBox1PropertyChange(PropertyChangeEvent e) {
        // TODO add your code here
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        setVisible(false);
    }

    private void showNetwork(List<String> geneLoci) {
        //TODO Run all of this on separate thread
        GeneNetwork network = GeneNetwork.getFromCBIO(geneLoci);

        network.annotateAll(IGV.getInstance().getAllTracks(false));

        //TODO Filter over multiple things
        String filt_el = (String) comboBox1.getSelectedItem();
        if (GeneNetwork.attribute_map.containsKey(filt_el)) {
            float min = rangeSlider1.getLowValue();
            float max = rangeSlider1.getHighValue();
            network.filterNodesRange(filt_el, min / 100, max / 100);
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


    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        comboBox1 = new JComboBox();
        rangeSlider1 = new RangeSlider();
        buttonBar = new JPanel();
        keepIsolated = new JCheckBox();
        okButton = new JButton();
        cancelButton = new JButton();
        helpButton = new JButton();

        //======== this ========
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.X_AXIS));

                //---- comboBox1 ----
                comboBox1.setToolTipText("Property by which to filter");
                contentPanel.add(comboBox1);

                //---- rangeSlider1 ----
                rangeSlider1.setLowValue(0);
                rangeSlider1.setHighValue(100);
                rangeSlider1.setPaintTicks(true);
                contentPanel.add(rangeSlider1);
            }
            dialogPane.add(contentPanel, BorderLayout.CENTER);

            //======== buttonBar ========
            {
                buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                buttonBar.setLayout(new GridBagLayout());
                ((GridBagLayout) buttonBar.getLayout()).columnWidths = new int[]{0, 85, 85, 80};
                ((GridBagLayout) buttonBar.getLayout()).columnWeights = new double[]{1.0, 0.0, 0.0, 0.0};

                //---- keepIsolated ----
                keepIsolated.setText("Keep Isolated Genes");
                buttonBar.add(keepIsolated, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0,
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
                buttonBar.add(okButton, new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0,
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
                buttonBar.add(cancelButton, new GridBagConstraints(2, 1, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 0, 5), 0, 0));

                //---- helpButton ----
                helpButton.setText("Help");
                helpButton.setVisible(false);
                buttonBar.add(helpButton, new GridBagConstraints(3, 1, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 0, 0), 0, 0));
            }
            dialogPane.add(buttonBar, BorderLayout.SOUTH);
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
    private JComboBox comboBox1;
    private RangeSlider rangeSlider1;
    private JPanel buttonBar;
    private JCheckBox keepIsolated;
    private JButton okButton;
    private JButton cancelButton;
    private JButton helpButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
