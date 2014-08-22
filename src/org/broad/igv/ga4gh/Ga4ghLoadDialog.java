/*
 * Created by JFormDesigner on Thu Aug 14 13:25:36 EDT 2014
 */

package org.broad.igv.ga4gh;

import org.broad.igv.ui.IGV;
import org.broad.igv.util.Pair;
import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.awt.event.*;
import java.util.Arrays;
import javax.swing.*;
import javax.swing.border.*;

/**
 * @author James Robinson
 */
public class Ga4ghLoadDialog extends JDialog {

    java.util.List<Pair<String, String>> idNamePairs;
    String selectedId;

    public Ga4ghLoadDialog(Frame owner, java.util.List<Pair<String, String>> idNamePairs) {
        super(owner);
        initComponents();

        this.idNamePairs = idNamePairs;
        String [] names = new String[idNamePairs.size()];
        for(int i=0; i<names.length; i++) {
            names[i] = idNamePairs.get(i).getSecond();
        }
        readsetSelectComboBox.setModel(new DefaultComboBoxModel(names));

    }

    private void loadButtonActionPerformed(ActionEvent e) {

        int idx = this.readsetSelectComboBox.getSelectedIndex();
        Pair<String, String> idName = idNamePairs.get(idx);
        setVisible(false);
        loadTrack(idName.getFirst(), idName.getSecond());
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        selectedId = null;
        setVisible(false);
    }


    private void loadTrack(String readsetId, String name) {

        ResourceLocator locator = new ResourceLocator(readsetId);
        locator.setName(name);
        locator.setType(Ga4ghAPIHelper.RESOURCE_TYPE);
        IGV.getInstance().loadTracks(Arrays.asList(locator));

    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner Evaluation license - James Robinson
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        readsetSelectComboBox = new JComboBox();
        label1 = new JLabel();
        buttonBar = new JPanel();
        loadButton = new JButton();
        cancelButton = new JButton();

        //======== this ========
        setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("GA4GH Prototype");
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));

            // JFormDesigner evaluation mark
           /* dialogPane.setBorder(new javax.swing.border.CompoundBorder(
                new javax.swing.border.TitledBorder(new javax.swing.border.EmptyBorder(0, 0, 0, 0),
                    "JFormDesigner Evaluation", javax.swing.border.TitledBorder.CENTER,
                    javax.swing.border.TitledBorder.BOTTOM, new java.awt.Font("Dialog", java.awt.Font.BOLD, 12),
                    java.awt.Color.red), dialogPane.getBorder())); dialogPane.addPropertyChangeListener(new java.beans.PropertyChangeListener(){public void propertyChange(java.beans.PropertyChangeEvent e){if("border".equals(e.getPropertyName()))throw new RuntimeException();}});
*/
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(null);
                contentPanel.add(readsetSelectComboBox);
                readsetSelectComboBox.setBounds(15, 60, 345, readsetSelectComboBox.getPreferredSize().height);

                //---- label1 ----
                label1.setText("Select a readset to load");
                contentPanel.add(label1);
                label1.setBounds(new Rectangle(new Point(15, 20), label1.getPreferredSize()));
            }
            dialogPane.add(contentPanel, BorderLayout.CENTER);

            //======== buttonBar ========
            {
                buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                buttonBar.setLayout(new GridBagLayout());
                ((GridBagLayout)buttonBar.getLayout()).columnWidths = new int[] {0, 85, 80};
                ((GridBagLayout)buttonBar.getLayout()).columnWeights = new double[] {1.0, 0.0, 0.0};

                //---- loadButton ----
                loadButton.setText("Load");
                loadButton.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        loadButtonActionPerformed(e);
                    }
                });
                buttonBar.add(loadButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
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
                buttonBar.add(cancelButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 0, 0), 0, 0));
            }
            dialogPane.add(buttonBar, BorderLayout.SOUTH);
        }
        contentPane.add(dialogPane, BorderLayout.CENTER);
        setSize(400, 225);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner Evaluation license - James Robinson
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JComboBox readsetSelectComboBox;
    private JLabel label1;
    private JPanel buttonBar;
    private JButton loadButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
