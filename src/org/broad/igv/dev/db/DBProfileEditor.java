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
 * Created by JFormDesigner on Thu Jul 05 13:43:44 EDT 2012
 */

package org.broad.igv.dev.db;

import org.broad.igv.feature.tribble.CodecFactory;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;

/**
 * @author User #2
 */
public class DBProfileEditor extends JDialog {
    public DBProfileEditor(Frame owner) {
        super(owner);
        initComponents();
        postInit();
    }

    public DBProfileEditor(Dialog owner) {
        super(owner);
        initComponents();
        postInit();
    }

    private void postInit() {
        //TODO Remember to add "." before extension when calling CodecFactory.getCodec
        DefaultComboBoxModel model = new DefaultComboBoxModel(CodecFactory.validExtensions.toArray(new String[0]));
        dataTypes.setModel(model);
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        panel4 = new JPanel();
        label4 = new JLabel();
        DBPath2 = new JTextField();
        panel1 = new JPanel();
        label1 = new JLabel();
        DBHost = new JTextField();
        panel6 = new JPanel();
        label6 = new JLabel();
        DBHost3 = new JTextField();
        panel5 = new JPanel();
        label5 = new JLabel();
        DBHost2 = new JTextField();
        panel2 = new JPanel();
        label2 = new JLabel();
        username = new JTextField();
        panel3 = new JPanel();
        label3 = new JLabel();
        password = new JPasswordField();
        checkBox1 = new JCheckBox();
        panel7 = new JPanel();
        label7 = new JLabel();
        DBPath3 = new JTextField();
        panel8 = new JPanel();
        label8 = new JLabel();
        dataTypes = new JComboBox();
        buttonBar = new JPanel();
        okButton = new JButton();
        cancelButton = new JButton();

        //======== this ========
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.Y_AXIS));

                //======== panel4 ========
                {
                    panel4.setLayout(new BoxLayout(panel4, BoxLayout.X_AXIS));

                    //---- label4 ----
                    label4.setText("Name:");
                    label4.setHorizontalAlignment(SwingConstants.LEFT);
                    label4.setMaximumSize(new Dimension(80, 16));
                    panel4.add(label4);

                    //---- DBPath2 ----
                    DBPath2.setMaximumSize(new Dimension(250, 28));
                    DBPath2.setPreferredSize(new Dimension(250, 28));
                    DBPath2.setToolTipText("mysql://my.awesomedb.com:8080/mytable");
                    DBPath2.setText("My Awesome DB");
                    panel4.add(DBPath2);
                }
                contentPanel.add(panel4);

                //======== panel1 ========
                {
                    panel1.setLayout(new BoxLayout(panel1, BoxLayout.X_AXIS));

                    //---- label1 ----
                    label1.setText("Host:");
                    label1.setHorizontalAlignment(SwingConstants.LEFT);
                    label1.setMaximumSize(new Dimension(80, 16));
                    panel1.add(label1);

                    //---- DBHost ----
                    DBHost.setMaximumSize(new Dimension(250, 28));
                    DBHost.setPreferredSize(new Dimension(250, 28));
                    DBHost.setToolTipText("my.awesomedb.com");
                    DBHost.setText("mysql://my.awesomedb.com");
                    panel1.add(DBHost);
                }
                contentPanel.add(panel1);

                //======== panel6 ========
                {
                    panel6.setLayout(new BoxLayout(panel6, BoxLayout.X_AXIS));

                    //---- label6 ----
                    label6.setText("Path:");
                    label6.setHorizontalAlignment(SwingConstants.LEFT);
                    label6.setMaximumSize(new Dimension(80, 16));
                    panel6.add(label6);

                    //---- DBHost3 ----
                    DBHost3.setMaximumSize(new Dimension(250, 28));
                    DBHost3.setPreferredSize(new Dimension(250, 28));
                    DBHost3.setToolTipText("my.awesomedb.com");
                    DBHost3.setText("hg19");
                    panel6.add(DBHost3);
                }
                contentPanel.add(panel6);

                //======== panel5 ========
                {
                    panel5.setLayout(new BoxLayout(panel5, BoxLayout.X_AXIS));

                    //---- label5 ----
                    label5.setText("Port (optional):");
                    label5.setHorizontalAlignment(SwingConstants.LEFT);
                    label5.setMaximumSize(new Dimension(80, 16));
                    panel5.add(label5);

                    //---- DBHost2 ----
                    DBHost2.setMaximumSize(new Dimension(250, 28));
                    DBHost2.setPreferredSize(new Dimension(250, 28));
                    DBHost2.setToolTipText("my.awesomedb.com");
                    DBHost2.setText("80");
                    panel5.add(DBHost2);
                }
                contentPanel.add(panel5);

                //======== panel2 ========
                {
                    panel2.setLayout(new BoxLayout(panel2, BoxLayout.X_AXIS));

                    //---- label2 ----
                    label2.setText("Username:");
                    label2.setHorizontalAlignment(SwingConstants.LEFT);
                    label2.setMaximumSize(new Dimension(80, 16));
                    label2.setMinimumSize(new Dimension(66, 16));
                    label2.setPreferredSize(new Dimension(80, 16));
                    panel2.add(label2);

                    //---- username ----
                    username.setMaximumSize(new Dimension(150, 28));
                    username.setPreferredSize(new Dimension(150, 28));
                    panel2.add(username);
                }
                contentPanel.add(panel2);

                //======== panel3 ========
                {
                    panel3.setLayout(new BoxLayout(panel3, BoxLayout.X_AXIS));

                    //---- label3 ----
                    label3.setText("Password (opt.):");
                    label3.setHorizontalAlignment(SwingConstants.LEFT);
                    label3.setMaximumSize(new Dimension(80, 16));
                    label3.setMinimumSize(new Dimension(66, 16));
                    label3.setPreferredSize(new Dimension(80, 16));
                    panel3.add(label3);

                    //---- password ----
                    password.setMaximumSize(new Dimension(150, 28));
                    password.setPreferredSize(new Dimension(150, 28));
                    password.setEnabled(false);
                    panel3.add(password);

                    //---- checkBox1 ----
                    checkBox1.setText("Save");
                    panel3.add(checkBox1);
                }
                contentPanel.add(panel3);

                //======== panel7 ========
                {
                    panel7.setLayout(new BoxLayout(panel7, BoxLayout.X_AXIS));

                    //---- label7 ----
                    label7.setText("Table Name:");
                    label7.setHorizontalAlignment(SwingConstants.LEFT);
                    label7.setMaximumSize(new Dimension(80, 16));
                    panel7.add(label7);

                    //---- DBPath3 ----
                    DBPath3.setMaximumSize(new Dimension(250, 28));
                    DBPath3.setPreferredSize(new Dimension(250, 28));
                    DBPath3.setToolTipText("mysql://my.awesomedb.com:8080/mytable");
                    DBPath3.setText("DataTable");
                    panel7.add(DBPath3);
                }
                contentPanel.add(panel7);

                //======== panel8 ========
                {
                    panel8.setLayout(new BoxLayout(panel8, BoxLayout.X_AXIS));

                    //---- label8 ----
                    label8.setText("Data Type:");
                    label8.setHorizontalAlignment(SwingConstants.LEFT);
                    label8.setMaximumSize(new Dimension(80, 16));
                    panel8.add(label8);

                    //---- dataTypes ----
                    dataTypes.setMaximumSize(new Dimension(250, 28));
                    dataTypes.setMinimumSize(new Dimension(96, 28));
                    panel8.add(dataTypes);
                }
                contentPanel.add(panel8);
            }
            dialogPane.add(contentPanel, BorderLayout.CENTER);

            //======== buttonBar ========
            {
                buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                buttonBar.setLayout(new GridBagLayout());
                ((GridBagLayout) buttonBar.getLayout()).columnWidths = new int[]{0, 85, 80};
                ((GridBagLayout) buttonBar.getLayout()).columnWeights = new double[]{1.0, 0.0, 0.0};

                //---- okButton ----
                okButton.setText("OK");
                buttonBar.add(okButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 0, 5), 0, 0));

                //---- cancelButton ----
                cancelButton.setText("Cancel");
                buttonBar.add(cancelButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0,
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
    private JPanel panel4;
    private JLabel label4;
    private JTextField DBPath2;
    private JPanel panel1;
    private JLabel label1;
    private JTextField DBHost;
    private JPanel panel6;
    private JLabel label6;
    private JTextField DBHost3;
    private JPanel panel5;
    private JLabel label5;
    private JTextField DBHost2;
    private JPanel panel2;
    private JLabel label2;
    private JTextField username;
    private JPanel panel3;
    private JLabel label3;
    private JPasswordField password;
    private JCheckBox checkBox1;
    private JPanel panel7;
    private JLabel label7;
    private JTextField DBPath3;
    private JPanel panel8;
    private JLabel label8;
    private JComboBox dataTypes;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
