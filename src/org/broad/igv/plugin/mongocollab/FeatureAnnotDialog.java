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
 * Created by JFormDesigner on Mon Nov 26 16:20:15 EST 2012
 */

package org.broad.igv.plugin.mongocollab;

import com.mongodb.DBCollection;
import org.apache.log4j.Logger;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.tribble.Feature;

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

    private String userName;
    private DBCollection collection;
    private DBFeature featDBObject;

    FeatureAnnotDialog(Frame owner, DBCollection collection, Feature feature) {
        super(owner);
        initComponents();

        this.userName = System.getProperty("user.name", "unknown");

        if(collection == null) throw new IllegalArgumentException("DBCollection must not be null");
        this.collection = collection;
        if(feature instanceof DBFeature){
            this.featDBObject = (DBFeature) feature;
        }else{
            this.featDBObject = DBFeature.create(feature);
        }

        initComponentData(this.featDBObject);
    }

    private void initComponentData(DBFeature feature) {

        setTitle("Add/Save Feature in " + this.collection.getFullName());

        userNameField.setText(userName);

        chrField.setText(feature.getChr());
        startField.setText("" + (feature.getStart() + 1));
        stopField.setText("" + feature.getEnd());

        scoreField.setText("" + feature.getScore());
        descField.setText("" + feature.getDescription());

        validate();
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        setVisible(false);
    }

    private DBFeature createDBObjectFromFields(){

        int start, end;
        float score;
        try{
            //Change from 1-based (user entered) to 0-based
            start = Integer.parseInt(startField.getText()) - 1;
            end = Integer.parseInt(stopField.getText());
            score = Float.parseFloat(scoreField.getText());
        }catch(NumberFormatException e){
            MessageUtils.showErrorMessage(e.getMessage(), e);
            return null;
        }

        if(start >= end){
            String msg = String.format("Start (%s) must be less than End (%s), but it's greater or equal", start, end);
            MessageUtils.showMessage(msg);
            return null;
        }

        this.featDBObject.setChr(chrField.getText());
        this.featDBObject.setStart(start);
        this.featDBObject.setEnd(end);
        this.featDBObject.setScore(score);
        this.featDBObject.setDescription(descField.getText());
        return this.featDBObject;
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
            //Find the track showing results, clear it to force a refresh
            for (FeatureTrack ft : IGV.getInstance().getFeatureTracks()) {
                //TODO This is a somewhat fragile way of identifying the corresponding track
                if(ft.getId() != null && ft.getId().equals(this.collection.getFullName())){
                    ft.clearPackedFeatures();
                }
            }
            IGV.getInstance().repaintDataPanels();
        }
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        panel1 = new JPanel();
        label4 = new JLabel();
        userNameField = new JTextField();
        panel4 = new JPanel();
        label7 = new JLabel();
        chrField = new JTextField();
        hSpacer3 = new JPanel(null);
        panel5 = new JPanel();
        label8 = new JLabel();
        startField = new JTextField();
        hSpacer4 = new JPanel(null);
        panel6 = new JPanel();
        label9 = new JLabel();
        stopField = new JTextField();
        hSpacer5 = new JPanel(null);
        panel7 = new JPanel();
        label10 = new JLabel();
        scoreField = new JTextField();
        hSpacer6 = new JPanel(null);
        panel8 = new JPanel();
        label11 = new JLabel();
        descField = new JTextField();
        buttonBar = new JPanel();
        okButton = new JButton();
        cancelButton = new JButton();

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

                //======== panel1 ========
                {
                    panel1.setLayout(new BoxLayout(panel1, BoxLayout.PAGE_AXIS));

                    //---- label4 ----
                    label4.setText("User");
                    label4.setHorizontalAlignment(SwingConstants.CENTER);
                    label4.setMaximumSize(new Dimension(80, 16));
                    label4.setAlignmentX(0.5F);
                    panel1.add(label4);

                    //---- userNameField ----
                    userNameField.setMaximumSize(new Dimension(200, 1000));
                    userNameField.setMinimumSize(new Dimension(100, 28));
                    userNameField.setPreferredSize(new Dimension(100, 28));
                    panel1.add(userNameField);
                }
                contentPanel.add(panel1);

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

                //---- hSpacer3 ----
                hSpacer3.setMinimumSize(new Dimension(20, 12));
                hSpacer3.setPreferredSize(new Dimension(20, 10));
                contentPanel.add(hSpacer3);

                //======== panel5 ========
                {
                    panel5.setMaximumSize(new Dimension(500, 1000));
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

                //---- hSpacer4 ----
                hSpacer4.setMinimumSize(new Dimension(20, 12));
                hSpacer4.setPreferredSize(new Dimension(20, 10));
                contentPanel.add(hSpacer4);

                //======== panel6 ========
                {
                    panel6.setMaximumSize(new Dimension(500, 1000));
                    panel6.setMinimumSize(new Dimension(100, 32));
                    panel6.setPreferredSize(new Dimension(100, 32));
                    panel6.setLayout(new BoxLayout(panel6, BoxLayout.Y_AXIS));

                    //---- label9 ----
                    label9.setText("End");
                    label9.setHorizontalAlignment(SwingConstants.LEFT);
                    label9.setMaximumSize(new Dimension(100, 16));
                    panel6.add(label9);
                    panel6.add(stopField);
                }
                contentPanel.add(panel6);

                //---- hSpacer5 ----
                hSpacer5.setMinimumSize(new Dimension(20, 12));
                hSpacer5.setPreferredSize(new Dimension(20, 10));
                contentPanel.add(hSpacer5);

                //======== panel7 ========
                {
                    panel7.setMaximumSize(new Dimension(46, 1000));
                    panel7.setMinimumSize(new Dimension(29, 44));
                    panel7.setPreferredSize(new Dimension(40, 44));
                    panel7.setLayout(new BoxLayout(panel7, BoxLayout.Y_AXIS));

                    //---- label10 ----
                    label10.setText("Score");
                    label10.setHorizontalAlignment(SwingConstants.LEFT);
                    label10.setMaximumSize(new Dimension(100, 16));
                    panel7.add(label10);
                    panel7.add(scoreField);
                }
                contentPanel.add(panel7);

                //---- hSpacer6 ----
                hSpacer6.setMinimumSize(new Dimension(20, 12));
                hSpacer6.setPreferredSize(new Dimension(20, 10));
                contentPanel.add(hSpacer6);

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
                    panel8.add(label11);
                    panel8.add(descField);
                }
                contentPanel.add(panel8);
            }
            dialogPane.add(contentPanel, BorderLayout.NORTH);

            //======== buttonBar ========
            {
                buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                buttonBar.setLayout(new GridBagLayout());
                ((GridBagLayout)buttonBar.getLayout()).columnWidths = new int[] {0, 85, 80};
                ((GridBagLayout)buttonBar.getLayout()).columnWeights = new double[] {1.0, 0.0, 0.0};

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
    private JPanel panel1;
    private JLabel label4;
    private JTextField userNameField;
    private JPanel panel4;
    private JLabel label7;
    private JTextField chrField;
    private JPanel hSpacer3;
    private JPanel panel5;
    private JLabel label8;
    private JTextField startField;
    private JPanel hSpacer4;
    private JPanel panel6;
    private JLabel label9;
    private JTextField stopField;
    private JPanel hSpacer5;
    private JPanel panel7;
    private JLabel label10;
    private JTextField scoreField;
    private JPanel hSpacer6;
    private JPanel panel8;
    private JLabel label11;
    private JTextField descField;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
