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
import org.broad.igv.util.Utilities;
import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.xml.parsers.ParserConfigurationException;
import java.awt.*;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;

/**
 * @author User #2
 */
public class DBProfileEditor extends JDialog {

    public DBProfileEditor(Frame owner, String initProfilePath) {
        super(owner);
        initComponents();
        postInit(initProfilePath);
    }

    public DBProfileEditor(Dialog owner, String initProfilePath) {
        super(owner);
        initComponents();
        postInit(initProfilePath);
    }

    private void postInit(String initProfilePath) {
        //TODO Remember to add "." before extension when calling CodecFactory.getCodec
        DefaultComboBoxModel model = new DefaultComboBoxModel(CodecFactory.validExtensions.toArray(new String[0]));
        dataType.setModel(model);

        if (initProfilePath != null) {
            InputStream profileStream;
            try {
                profileStream = new FileInputStream(initProfilePath);
                Document document = Utilities.createDOMDocumentFromXmlStream(profileStream);
                NamedNodeMap dbAttrs = document.getAttributes();

                String fullHost = "";
                String subprotocol = Utilities.getNullSafe(dbAttrs, "subprotocol");
                String host = Utilities.getNullSafe(dbAttrs, "host");
                if (subprotocol != null && host != null) {
                    fullHost = subprotocol + "://" + host;
                }
                DBPath.setText(fullHost);
                port.setText(Utilities.getNullSafe(dbAttrs, "port"));
                username.setText(Utilities.getNullSafe(dbAttrs, "username"));
                password.setText(Utilities.getNullSafe(dbAttrs, "password"));

                //TODO Can have more than 1 table, for now just take first
                NodeList tables = document.getElementsByTagName("table");
                Node table = tables.item(0);
                NamedNodeMap tableAttrs = table.getAttributes();
                tableName.setText(Utilities.getNullSafe(tableAttrs, "name"));
                chromField.setText(Utilities.getNullSafe(tableAttrs, "chromoColName"));

                posStartField.setText(Utilities.getNullSafe(tableAttrs, "posStartColName"));
                posEndField.setText(Utilities.getNullSafe(tableAttrs, "posEndColName"));

                startColField.setText(Utilities.getNullSafe(tableAttrs, "startColIndex"));
                endColField.setText(Utilities.getNullSafe(tableAttrs, "endColIndex"));
                binColField.setText(Utilities.getNullSafe(tableAttrs, "binColName"));

            } catch (FileNotFoundException e) {
                e.printStackTrace();
            } catch (ParserConfigurationException e) {
                e.printStackTrace();
            } catch (SAXException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            }

        }
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        panel4 = new JPanel();
        label4 = new JLabel();
        nameComboBox = new JComboBox();
        panel1 = new JPanel();
        label1 = new JLabel();
        DBHost = new JTextField();
        panel6 = new JPanel();
        label6 = new JLabel();
        DBPath = new JTextField();
        panel5 = new JPanel();
        label5 = new JLabel();
        port = new JTextField();
        panel2 = new JPanel();
        label2 = new JLabel();
        username = new JTextField();
        panel3 = new JPanel();
        label3 = new JLabel();
        password = new JPasswordField();
        checkBox1 = new JCheckBox();
        panel7 = new JPanel();
        label7 = new JLabel();
        tableName = new JTextField();
        panel8 = new JPanel();
        label8 = new JLabel();
        dataType = new JComboBox();
        separator1 = new JSeparator();
        panel9 = new JPanel();
        label9 = new JLabel();
        chromField = new JTextField();
        panel10 = new JPanel();
        label10 = new JLabel();
        posStartField = new JTextField();
        panel13 = new JPanel();
        label13 = new JLabel();
        posEndField = new JTextField();
        panel11 = new JPanel();
        label11 = new JLabel();
        startColField = new JTextField();
        panel12 = new JPanel();
        label12 = new JLabel();
        endColField = new JTextField();
        panel14 = new JPanel();
        label14 = new JLabel();
        binColField = new JTextField();
        buttonBar = new JPanel();
        saveButton = new JButton();
        okButton = new JButton();
        cancelButton = new JButton();

        //======== this ========
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setPreferredSize(new Dimension(367, 500));
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

                    //---- nameComboBox ----
                    nameComboBox.setEditable(true);
                    nameComboBox.setMaximumRowCount(100);
                    panel4.add(nameComboBox);
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
                    DBHost.setText("mysql://my.db.com");
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

                    //---- DBPath ----
                    DBPath.setMaximumSize(new Dimension(250, 28));
                    DBPath.setPreferredSize(new Dimension(250, 28));
                    DBPath.setToolTipText("my.awesomedb.com");
                    DBPath.setText("hg19");
                    panel6.add(DBPath);
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

                    //---- port ----
                    port.setMaximumSize(new Dimension(250, 28));
                    port.setPreferredSize(new Dimension(250, 28));
                    port.setToolTipText("my.awesomedb.com");
                    port.setText("80");
                    panel5.add(port);
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
                    panel3.setMaximumSize(new Dimension(400, 28));
                    panel3.setLayout(new BoxLayout(panel3, BoxLayout.X_AXIS));

                    //---- label3 ----
                    label3.setText("Password (opt.):");
                    label3.setMaximumSize(new Dimension(120, 16));
                    label3.setMinimumSize(new Dimension(66, 16));
                    label3.setPreferredSize(new Dimension(80, 16));
                    label3.setHorizontalAlignment(SwingConstants.LEFT);
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

                    //---- tableName ----
                    tableName.setMaximumSize(new Dimension(250, 28));
                    tableName.setPreferredSize(new Dimension(250, 28));
                    tableName.setToolTipText("mysql://my.awesomedb.com:8080/mytable");
                    tableName.setText("DataTable");
                    panel7.add(tableName);
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

                    //---- dataType ----
                    dataType.setMaximumSize(new Dimension(250, 28));
                    dataType.setMinimumSize(new Dimension(96, 28));
                    panel8.add(dataType);
                }
                contentPanel.add(panel8);

                //---- separator1 ----
                separator1.setPreferredSize(new Dimension(0, 1));
                contentPanel.add(separator1);

                //======== panel9 ========
                {
                    panel9.setLayout(new BoxLayout(panel9, BoxLayout.X_AXIS));

                    //---- label9 ----
                    label9.setText("Chromosome Column Name:");
                    label9.setHorizontalAlignment(SwingConstants.LEFT);
                    label9.setMaximumSize(new Dimension(80, 16));
                    panel9.add(label9);

                    //---- chromField ----
                    chromField.setMaximumSize(new Dimension(250, 28));
                    chromField.setPreferredSize(new Dimension(250, 28));
                    chromField.setText("chrom");
                    panel9.add(chromField);
                }
                contentPanel.add(panel9);

                //======== panel10 ========
                {
                    panel10.setLayout(new BoxLayout(panel10, BoxLayout.X_AXIS));

                    //---- label10 ----
                    label10.setText("Position Start Column Name:");
                    label10.setHorizontalAlignment(SwingConstants.LEFT);
                    label10.setMaximumSize(new Dimension(80, 16));
                    panel10.add(label10);

                    //---- posStartField ----
                    posStartField.setMaximumSize(new Dimension(250, 28));
                    posStartField.setPreferredSize(new Dimension(250, 28));
                    posStartField.setText("txStart");
                    panel10.add(posStartField);
                }
                contentPanel.add(panel10);

                //======== panel13 ========
                {
                    panel13.setLayout(new BoxLayout(panel13, BoxLayout.X_AXIS));

                    //---- label13 ----
                    label13.setText("Position End Column Name:");
                    label13.setHorizontalAlignment(SwingConstants.LEFT);
                    label13.setMaximumSize(new Dimension(80, 16));
                    panel13.add(label13);

                    //---- posEndField ----
                    posEndField.setMaximumSize(new Dimension(250, 28));
                    posEndField.setPreferredSize(new Dimension(250, 28));
                    posEndField.setText("txEnd");
                    panel13.add(posEndField);
                }
                contentPanel.add(panel13);

                //======== panel11 ========
                {
                    panel11.setLayout(new BoxLayout(panel11, BoxLayout.X_AXIS));

                    //---- label11 ----
                    label11.setText("Data start column index (opt.):");
                    label11.setHorizontalAlignment(SwingConstants.LEFT);
                    label11.setMaximumSize(new Dimension(80, 16));
                    panel11.add(label11);

                    //---- startColField ----
                    startColField.setMaximumSize(new Dimension(250, 28));
                    startColField.setPreferredSize(new Dimension(250, 28));
                    startColField.setToolTipText("Starting column index from which to read data (1-based)");
                    startColField.setText("1");
                    panel11.add(startColField);
                }
                contentPanel.add(panel11);

                //======== panel12 ========
                {
                    panel12.setLayout(new BoxLayout(panel12, BoxLayout.X_AXIS));

                    //---- label12 ----
                    label12.setText("Data end column index (opt.):");
                    label12.setHorizontalAlignment(SwingConstants.LEFT);
                    label12.setMaximumSize(new Dimension(80, 16));
                    panel12.add(label12);

                    //---- endColField ----
                    endColField.setMaximumSize(new Dimension(250, 28));
                    endColField.setPreferredSize(new Dimension(250, 28));
                    endColField.setToolTipText("Last column (1-based, inclusive end) from which to read data");
                    panel12.add(endColField);
                }
                contentPanel.add(panel12);

                //======== panel14 ========
                {
                    panel14.setLayout(new BoxLayout(panel14, BoxLayout.X_AXIS));

                    //---- label14 ----
                    label14.setText("Bin column Name (opt.):");
                    label14.setHorizontalAlignment(SwingConstants.LEFT);
                    label14.setMaximumSize(new Dimension(80, 16));
                    panel14.add(label14);

                    //---- binColField ----
                    binColField.setMaximumSize(new Dimension(250, 28));
                    binColField.setPreferredSize(new Dimension(250, 28));
                    binColField.setToolTipText("Some databases use binning to speed up querying");
                    panel14.add(binColField);
                }
                contentPanel.add(panel14);
            }
            dialogPane.add(contentPanel, BorderLayout.CENTER);

            //======== buttonBar ========
            {
                buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                buttonBar.setLayout(new GridBagLayout());
                ((GridBagLayout) buttonBar.getLayout()).columnWidths = new int[]{0, 85, 80};
                ((GridBagLayout) buttonBar.getLayout()).columnWeights = new double[]{1.0, 0.0, 0.0};

                //---- saveButton ----
                saveButton.setText("Save Profile");
                buttonBar.add(saveButton, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 0, 5), 0, 0));

                //---- okButton ----
                okButton.setText("Load Data");
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
    private JComboBox nameComboBox;
    private JPanel panel1;
    private JLabel label1;
    private JTextField DBHost;
    private JPanel panel6;
    private JLabel label6;
    private JTextField DBPath;
    private JPanel panel5;
    private JLabel label5;
    private JTextField port;
    private JPanel panel2;
    private JLabel label2;
    private JTextField username;
    private JPanel panel3;
    private JLabel label3;
    private JPasswordField password;
    private JCheckBox checkBox1;
    private JPanel panel7;
    private JLabel label7;
    private JTextField tableName;
    private JPanel panel8;
    private JLabel label8;
    private JComboBox dataType;
    private JSeparator separator1;
    private JPanel panel9;
    private JLabel label9;
    private JTextField chromField;
    private JPanel panel10;
    private JLabel label10;
    private JTextField posStartField;
    private JPanel panel13;
    private JLabel label13;
    private JTextField posEndField;
    private JPanel panel11;
    private JLabel label11;
    private JTextField startColField;
    private JPanel panel12;
    private JLabel label12;
    private JTextField endColField;
    private JPanel panel14;
    private JLabel label14;
    private JTextField binColField;
    private JPanel buttonBar;
    private JButton saveButton;
    private JButton okButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
