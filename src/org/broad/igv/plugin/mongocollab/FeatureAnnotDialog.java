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
 * Created by JFormDesigner on Mon Nov 26 16:20:15 EST 2012
 */

package org.broad.igv.plugin.mongocollab;

import com.mongodb.DBCollection;
import org.apache.log4j.Logger;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.color.ColorChooserPanel;
import org.broad.igv.ui.util.MessageUtils;
import htsjdk.tribble.Feature;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;

/**
 * Dialog for inserting a feature annotation
 * into collaborative database
 */
public class FeatureAnnotDialog extends JDialog {

    private static Logger log = Logger.getLogger(MongoCollabPlugin.class);

    private DBCollection collection;
    private DBFeature dbFeat;

    FeatureAnnotDialog(Frame owner, DBCollection collection, Feature feature) {
        super(owner);
        initComponents();

        if(collection == null) throw new IllegalArgumentException("DBCollection must not be null");
        this.collection = collection;
        if(feature instanceof DBFeature){
            this.dbFeat = (DBFeature) feature;
        }else{
            this.dbFeat = DBFeature.create(feature);
        }
        initComponentData(this.dbFeat);
    }

    private void initComponentData(DBFeature feature) {

        setTitle("Add/Save Feature in " + this.collection.getFullName());

        chrField.setText(feature.getChr());
        startField.setText("" + (feature.getStart() + 1));
        stopField.setText("" + feature.getEnd());

        nameField.setText(feature.getName());
        descField.setText("" + feature.getDescription());

        colorField.setSelectedColor(feature.getColor());

        if(this.dbFeat.get_id() == null){
            delButton.setVisible(false);
        }

        validate();
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        setVisible(false);
    }

    private DBFeature createDBObjectFromFields(){

        int start, end;
        try{
            //Change from 1-based (user entered) to 0-based
            start = Integer.parseInt(startField.getText()) - 1;
            end = Integer.parseInt(stopField.getText());
        }catch(NumberFormatException e){
            MessageUtils.showErrorMessage(e.getMessage(), e);
            return null;
        }

        if(start >= end){
            String msg = String.format("Start (%s) must be less than End (%s), but it's greater or equal", start, end);
            MessageUtils.showMessage(msg);
            return null;
        }

        this.dbFeat.setChr(chrField.getText());
        this.dbFeat.setStart(start);
        this.dbFeat.setEnd(end);
        this.dbFeat.setColor(colorField.getSelectedColor());
        this.dbFeat.setDescription(descField.getText());
        this.dbFeat.setName(nameField.getText());

        return this.dbFeat;
    }

    /**
     * Save information to MongoDB
     *
     * @param e
     */
    private void okButtonActionPerformed(ActionEvent e) {

        DBFeature featDBObject = createDBObjectFromFields();
        String errorMessage = MongoCollabPlugin.saveFeature(collection, featDBObject);

        if (errorMessage != null) {
            MessageUtils.showErrorMessage(errorMessage, new IOException(errorMessage));
        } else {
            setVisible(false);
            refreshTrack();
        }

    }

    private void refreshTrack(){
        //Find the track showing results, clear it to force a refresh
        for (FeatureTrack ft : IGV.getInstance().getFeatureTracks()) {
            //TODO This is a somewhat fragile way of identifying the corresponding track
            if(ft.getId() != null && ft.getId().equals(this.collection.getFullName())){
                ft.clearPackedFeatures();
            }
        }
        IGV.getInstance().repaintDataPanels();
    }

    private void delButtonActionPerformed(ActionEvent e) {
        MongoCollabPlugin.removeFeature(collection, this.dbFeat);
        setVisible(false);
        refreshTrack();
    }

    private void createUIComponents() {
        // TODO: add custom component creation code here
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        panel4 = new JPanel();
        label7 = new JLabel();
        chrField = new JTextField();
        panel5 = new JPanel();
        label8 = new JLabel();
        startField = new JTextField();
        panel6 = new JPanel();
        label9 = new JLabel();
        stopField = new JTextField();
        panel7 = new JPanel();
        color = new JLabel();
        colorField = new ColorChooserPanel();
        panel9 = new JPanel();
        name = new JLabel();
        nameField = new JTextField();
        panel8 = new JPanel();
        label11 = new JLabel();
        descField = new JTextField();
        buttonBar = new JPanel();
        okButton = new JButton();
        cancelButton = new JButton();
        delButton = new JButton();

        //======== this ========
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setPreferredSize(new Dimension(600, 150));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.X_AXIS));

                //======== panel4 ========
                {
                    panel4.setMaximumSize(new Dimension(46, 1000));
                    panel4.setLayout(new BoxLayout(panel4, BoxLayout.Y_AXIS));

                    //---- label7 ----
                    label7.setText("Chr");
                    label7.setHorizontalAlignment(SwingConstants.LEFT);
                    label7.setVerticalAlignment(SwingConstants.TOP);
                    label7.setLabelFor(chrField);
                    panel4.add(label7);
                    panel4.add(chrField);
                }
                contentPanel.add(panel4);

                //======== panel5 ========
                {
                    panel5.setMaximumSize(new Dimension(100, 1000));
                    panel5.setMinimumSize(new Dimension(80, 32));
                    panel5.setPreferredSize(new Dimension(80, 32));
                    panel5.setLayout(new BoxLayout(panel5, BoxLayout.Y_AXIS));

                    //---- label8 ----
                    label8.setText("Start");
                    label8.setHorizontalAlignment(SwingConstants.LEFT);
                    label8.setMaximumSize(new Dimension(100, 16));
                    panel5.add(label8);
                    panel5.add(startField);
                }
                contentPanel.add(panel5);

                //======== panel6 ========
                {
                    panel6.setMaximumSize(new Dimension(100, 1000));
                    panel6.setMinimumSize(new Dimension(80, 32));
                    panel6.setPreferredSize(new Dimension(80, 32));
                    panel6.setLayout(new BoxLayout(panel6, BoxLayout.Y_AXIS));

                    //---- label9 ----
                    label9.setText("End");
                    label9.setHorizontalAlignment(SwingConstants.LEFT);
                    label9.setMaximumSize(new Dimension(100, 16));
                    panel6.add(label9);
                    panel6.add(stopField);
                }
                contentPanel.add(panel6);

                //======== panel7 ========
                {
                    panel7.setMaximumSize(new Dimension(80, 1000));
                    panel7.setMinimumSize(new Dimension(60, 44));
                    panel7.setPreferredSize(new Dimension(60, 44));
                    panel7.setToolTipText("Click to edit");
                    panel7.setLayout(new BoxLayout(panel7, BoxLayout.Y_AXIS));

                    //---- color ----
                    color.setText("Color");
                    color.setHorizontalAlignment(SwingConstants.LEFT);
                    color.setMaximumSize(new Dimension(100, 16));
                    panel7.add(color);

                    //---- colorField ----
                    colorField.setSelectedColor(DBFeature.DEFAULT_COLOR);
                    panel7.add(colorField);
                }
                contentPanel.add(panel7);

                //======== panel9 ========
                {
                    panel9.setMaximumSize(new Dimension(200, 1000));
                    panel9.setMinimumSize(new Dimension(100, 44));
                    panel9.setPreferredSize(new Dimension(100, 44));
                    panel9.setLayout(new BoxLayout(panel9, BoxLayout.Y_AXIS));

                    //---- name ----
                    name.setText("Name");
                    name.setHorizontalAlignment(SwingConstants.LEFT);
                    name.setMaximumSize(new Dimension(100, 16));
                    panel9.add(name);
                    panel9.add(nameField);
                }
                contentPanel.add(panel9);

                //======== panel8 ========
                {
                    panel8.setMaximumSize(new Dimension(500, 1000));
                    panel8.setMinimumSize(new Dimension(100, 32));
                    panel8.setPreferredSize(new Dimension(100, 32));
                    panel8.setLayout(new BoxLayout(panel8, BoxLayout.Y_AXIS));

                    //---- label11 ----
                    label11.setText("Description");
                    label11.setHorizontalAlignment(SwingConstants.LEFT);
                    label11.setMaximumSize(new Dimension(100, 16));
                    label11.setLabelFor(descField);
                    label11.setAlignmentX(1.0F);
                    panel8.add(label11);

                    //---- descField ----
                    descField.setPreferredSize(new Dimension(100, 28));
                    descField.setMinimumSize(new Dimension(100, 28));
                    panel8.add(descField);
                }
                contentPanel.add(panel8);
            }
            dialogPane.add(contentPanel, BorderLayout.NORTH);

            //======== buttonBar ========
            {
                buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                buttonBar.setLayout(new GridBagLayout());
                ((GridBagLayout)buttonBar.getLayout()).columnWidths = new int[] {0, 85, 85, 0};
                ((GridBagLayout)buttonBar.getLayout()).columnWeights = new double[] {1.0, 0.0, 0.0, 0.0};

                //---- okButton ----
                okButton.setText("Save");
                okButton.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        okButtonActionPerformed(e);
                    }
                });
                buttonBar.add(okButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
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
                    new Insets(0, 0, 0, 5), 0, 0));

                //---- delButton ----
                delButton.setText("Delete");
                delButton.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        delButtonActionPerformed(e);
                    }
                });
                buttonBar.add(delButton, new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 0, 0), 0, 0));
            }
            dialogPane.add(buttonBar, BorderLayout.SOUTH);
        }
        contentPane.add(dialogPane, BorderLayout.CENTER);
        setSize(700, 160);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JPanel panel4;
    private JLabel label7;
    private JTextField chrField;
    private JPanel panel5;
    private JLabel label8;
    private JTextField startField;
    private JPanel panel6;
    private JLabel label9;
    private JTextField stopField;
    private JPanel panel7;
    private JLabel color;
    private ColorChooserPanel colorField;
    private JPanel panel9;
    private JLabel name;
    private JTextField nameField;
    private JPanel panel8;
    private JLabel label11;
    private JTextField descField;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    private JButton delButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
