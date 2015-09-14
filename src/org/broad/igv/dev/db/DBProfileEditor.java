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
 * Created by JFormDesigner on Thu Jul 05 13:43:44 EDT 2012
 */

package org.broad.igv.dev.db;

import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.Utilities;
import org.w3c.dom.Document;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.xml.bind.Marshaller;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Dialog for creating/editing a profile for reading data from a SQL database
 * TODO Bind beans the easy way instead of manually
 * @author jacob
 */
public class DBProfileEditor extends JDialog {

    private DBProfile profile = new DBProfile();
    private String profilePath = null;

    private Map<String, DBProfile.DBTable> allTableNames = new LinkedHashMap<String, DBProfile.DBTable>();

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
        tableFields = new JComponent[]{chromField, posStartField, posEndField, startColField, endColField, binColField, dataType};
        assert initProfilePath != null;

        //TODO Remember to add "." before extension when calling CodecFactory.getCodec
        DefaultComboBoxModel model = new DefaultComboBoxModel(CodecFactory.validExtensions.toArray(new String[0]));
        dataType.setModel(model);

        this.profilePath = initProfilePath;
        DBProfile.DBTable initTable = null;
        File initProfileFile = new File(initProfilePath);

        if (initProfileFile.exists()) {
            //Editing existing profile
            profile = DBProfile.parseProfile(initProfilePath);

            DBSubprotocol.setSelectedItem(profile.getSubprotocol());
            DBName.setText(profile.getName());
            DBHost.setText(profile.getHost());
            DBPath.setText(profile.getPath());

            port.setText(profile.getPort());
            username.setText(profile.getUsername());
            password.setText(profile.getPassword());

            initTableNameList();
            initTable = profile.getTableList().get(0);

        }else{
            //Creating new profile
            profile = new DBProfile();
            this.profilePath = initProfilePath;
            if (!this.profilePath.endsWith(".dbxml")) {
                this.profilePath += ".dbxml";
                initTable = null;
            }
        }


        populateTableFieldValues(initTable);

        attachTableFieldListeners();
    }

    /**
     * Fill in the combo box of table names
     */
    private void initTableNameList() {
        List<DBProfile.DBTable> tableList = profile.getTableList();
        for(DBProfile.DBTable table: tableList){
            allTableNames.put(table.getName(), table);
        }

        tableName.setModel(new VetoableComboBoxModel(allTableNames.keySet().toArray(new String[tableList.size()])));

//        tableName.addItemListener(new ItemListener() {
//            @Override
//            public void itemStateChanged(ItemEvent e) {
//                //If selected table not done, don't let user change
//                if(e.getStateChange() == ItemEvent.DESELECTED){
//                    String tabName = (String) e.getItem();
//                    if(!checkTableAndWarn(tabName)){
//                        tableName.setSelectedItem(tabName);
//                        populateTableFieldValues(allTableNames.get(tabName));
//                    }
//                }
//            }
//        });
    }

    /**
     * Save text inputs to {@link #profile}
     * Doesn't write to disk
     */
    private void saveDBInputs(){
        profile.setSubprotocol(DBSubprotocol.getSelectedItem().toString());
        profile.setName(DBName.getText());
        profile.setHost(DBHost.getText());
        profile.setPath(DBPath.getText());

        profile.setPort(port.getText());
        profile.setUsername(username.getText());

        char[] cpw = password.getPassword();
        if(cpw != null && cpw.length > 0){
            String spw = new String(cpw);
            profile.setPassword(spw);
            spw = null;
            Arrays.fill(cpw, (char) 0);
        }

    }

    private DBProfile.DBTable getSelectedTable(){
        String selectedTableName = (String) this.tableName.getSelectedItem();
        return allTableNames.get(selectedTableName);
    }

    /**
     * Save text inputs to the appropriate table
     * Doesn't write to disk
     * @param table
     * @return true if data saved
     */
    private boolean saveTableInput(DBProfile.DBTable table){

        if(table == null) return false;

        int startColIndex = -1;
        int endColIndex = -1;

        try{
            startColIndex = Integer.parseInt(startColField.getText());
            endColIndex = Integer.parseInt(endColField.getText());
        }catch (NumberFormatException e){
            MessageUtils.showErrorMessage("Entry must be a valid integer: " + e.getMessage(), e);
            return false;
        }

        table.setStartColIndex(startColIndex);
        table.setEndColIndex(endColIndex);

        table.setChromoColName(chromField.getText());
        table.setPosStartColName(posStartField.getText());
        table.setPosEndColName(posEndField.getText());
        table.setBinColName(binColField.getText());
        table.setFormat((String) dataType.getSelectedItem());

        return true;
    }

    /**
     * We save the text field inputs to the relevant {@link #profile#table}
     * whenever focus is lost.
     */
    private void attachTableFieldListeners(){

        FocusListener tableFocusListener = new FocusListener() {
            @Override
            public void focusGained(FocusEvent e) {
                //saveTableInput(getSelectedTable());
            }

            @Override
            public void focusLost(FocusEvent e) {
                saveTableInput(getSelectedTable());
            }
        };

        for(JComponent tableField: tableFields){
            tableField.addFocusListener(tableFocusListener);
        }

    }

    private void populateTableFieldValues(DBProfile.DBTable table){

        boolean enabled = table != null;
        for(JComponent tableField: tableFields){
            tableField.setEnabled(enabled);
        }

        if(table == null) return;


        chromField.setText(table.getChromoColName());

        posStartField.setText(table.getPosStartColName());
        posEndField.setText(table.getPosEndColName());

        startColField.setText("" + table.getStartColIndex());
        endColField.setText("" + table.getEndColIndex());
        binColField.setText(table.getBinColName());

        dataType.setSelectedItem(table.getFormat());
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        setVisible(false);
    }

    private void tableNameActionPerformed(ActionEvent e) {
        populateTableFieldValues(getSelectedTable());
    }

    private void saveButtonActionPerformed(ActionEvent e) {
        saveDBInputs();

        if(!checkDBInputs()){
            return;
        }

        //Save to disk
        FileWriter fileWriter = null;
        try {
            // Create a DOM document
            DocumentBuilder documentBuilder = DocumentBuilderFactory.newInstance().newDocumentBuilder();
            Document document = documentBuilder.newDocument();
            document.setStrictErrorChecking(false);

            Marshaller m = DBProfile.getJAXBContext().createMarshaller();
            m.marshal(profile, document);

            String xmlString = Utilities.getString(document);
            fileWriter = new FileWriter(this.profilePath);
            fileWriter.write(xmlString);

            this.setVisible(false);
        }catch(Exception ex){
            ex.printStackTrace();
        }finally {
            if (fileWriter != null) {
                try {
                    fileWriter.close();
                } catch (IOException e1) {
                    e1.printStackTrace();
                }
            }

        }


    }

    /**
     * Checks that all fields have been filled out,
     * and if not, warns the user
     * @return true if all required fields have been filled out, else false
     */
    private boolean checkDBInputs() {
        List<String> missingDBFields = this.profile.checkMissingValues();
        if(missingDBFields.size() > 0){
            String msg = String.format("Please fill in all required database fields");
            MessageUtils.showMessage(msg);
            return false;
        }

        List<DBProfile.DBTable> tables = this.profile.getTableList();
        if(tables.size() == 0){
            MessageUtils.showMessage("You must add at least one table");
            return false;
        }

        for(DBProfile.DBTable table: tables){
            if(!checkTableAndWarn(table.getName())){
                return false;
            }
        }
        return true;
    }

    private void addNewTableButtonActionPerformed(ActionEvent e) {
        saveDBInputs();

        if(!checkTableAndWarn((String) tableName.getSelectedItem())){
            return;
        }

        String strTableName = newTableNameField.getText();
        if(strTableName.length() == 0){
            MessageUtils.showMessage("Please enter a table name");
            return;
        }
        DBProfile.DBTable newTable = new DBProfile.DBTable(profile.getDBLocator(), strTableName);

        profile.addTable(newTable);
        initTableNameList();

        tableName.setSelectedItem(strTableName);
        newTableNameField.setText("");
    }

    /**
     * Verify that all fields for the table
     * have been filled out. If the specified {@code tableName} is null,
     * return true
     * @return
     */
    private boolean checkTable(String tableName) {
        if(tableName == null) return true;
        DBProfile.DBTable table = allTableNames.get(tableName);
        List<String> missingFields = table.checkMissingValues();
        return missingFields.size() == 0;
    }

    private boolean checkTableAndWarn(String tableName){
        if(!checkTable(tableName)){
            MessageUtils.showMessage("Please fill in all required table fields in table " + tableName);
            return false;
        }
        return true;
    }

    /**
     * ComboBoxModel which doesn't change if fields don't pass verification
     */
    private class VetoableComboBoxModel extends DefaultComboBoxModel{

        private VetoableComboBoxModel(Object[] objects){
            super(objects);
        }

        @Override
        public void setSelectedItem(Object anItem) {
            Object oldItem = getSelectedItem();
            if(!checkTableAndWarn((String) oldItem)){
                return;
            }
            super.setSelectedItem(anItem);

        }
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        panel16 = new JPanel();
        label16 = new JLabel();
        DBSubprotocol = new JComboBox();
        panel4 = new JPanel();
        label4 = new JLabel();
        DBName = new JTextField();
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
        separator1 = new JSeparator();
        panel15 = new JPanel();
        label15 = new JLabel();
        newTableNameField = new JTextField();
        addNewTableButton = new JButton();
        separator2 = new JSeparator();
        panel7 = new JPanel();
        label7 = new JLabel();
        tableName = new JComboBox();
        panel8 = new JPanel();
        label8 = new JLabel();
        dataType = new JComboBox();
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
            dialogPane.setPreferredSize(new Dimension(367, 550));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.Y_AXIS));

                //======== panel16 ========
                {
                    panel16.setLayout(new BoxLayout(panel16, BoxLayout.X_AXIS));

                    //---- label16 ----
                    label16.setText("Subprotocol:");
                    label16.setHorizontalAlignment(SwingConstants.LEFT);
                    label16.setMaximumSize(new Dimension(80, 16));
                    panel16.add(label16);

                    //---- DBSubprotocol ----
                    DBSubprotocol.setMaximumSize(new Dimension(250, 28));
                    DBSubprotocol.setPreferredSize(new Dimension(250, 28));
                    DBSubprotocol.setModel(new DefaultComboBoxModel(new String[] {
                        "mysql",
                        "sqlite"
                    }));
                    DBSubprotocol.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            tableNameActionPerformed(e);
                        }
                    });
                    panel16.add(DBSubprotocol);
                }
                contentPanel.add(panel16);

                //======== panel4 ========
                {
                    panel4.setMaximumSize(new Dimension(330, 28));
                    panel4.setLayout(new BoxLayout(panel4, BoxLayout.X_AXIS));

                    //---- label4 ----
                    label4.setText("Name:");
                    label4.setHorizontalAlignment(SwingConstants.LEFT);
                    label4.setMaximumSize(new Dimension(80, 16));
                    panel4.add(label4);

                    //---- DBName ----
                    DBName.setEditable(true);
                    panel4.add(DBName);
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
                    DBPath.setToolTipText("myfolder");
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

                //---- separator1 ----
                separator1.setPreferredSize(new Dimension(0, 1));
                contentPanel.add(separator1);

                //======== panel15 ========
                {
                    panel15.setLayout(new BoxLayout(panel15, BoxLayout.X_AXIS));

                    //---- label15 ----
                    label15.setText("New Table Name:");
                    panel15.add(label15);

                    //---- newTableNameField ----
                    newTableNameField.setToolTipText("Name for new table");
                    newTableNameField.setMaximumSize(new Dimension(250, 28));
                    panel15.add(newTableNameField);

                    //---- addNewTableButton ----
                    addNewTableButton.setText("Add");
                    addNewTableButton.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            addNewTableButtonActionPerformed(e);
                        }
                    });
                    panel15.add(addNewTableButton);
                }
                contentPanel.add(panel15);

                //---- separator2 ----
                separator2.setPreferredSize(new Dimension(0, 1));
                contentPanel.add(separator2);

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
                    tableName.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            tableNameActionPerformed(e);
                        }
                    });
                    panel7.add(tableName);
                }
                contentPanel.add(panel7);

                //======== panel8 ========
                {
                    panel8.setLayout(new BoxLayout(panel8, BoxLayout.X_AXIS));

                    //---- label8 ----
                    label8.setText("Format:");
                    label8.setHorizontalAlignment(SwingConstants.LEFT);
                    label8.setMaximumSize(new Dimension(80, 16));
                    panel8.add(label8);

                    //---- dataType ----
                    dataType.setMaximumSize(new Dimension(250, 28));
                    dataType.setMinimumSize(new Dimension(96, 28));
                    panel8.add(dataType);
                }
                contentPanel.add(panel8);

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
                    panel13.add(posEndField);
                }
                contentPanel.add(panel13);

                //======== panel11 ========
                {
                    panel11.setLayout(new BoxLayout(panel11, BoxLayout.X_AXIS));

                    //---- label11 ----
                    label11.setText("Data start column index:");
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
                    label12.setText("Data end column index:");
                    label12.setHorizontalAlignment(SwingConstants.LEFT);
                    label12.setMaximumSize(new Dimension(80, 16));
                    panel12.add(label12);

                    //---- endColField ----
                    endColField.setMaximumSize(new Dimension(250, 28));
                    endColField.setPreferredSize(new Dimension(250, 28));
                    endColField.setToolTipText("Last column (1-based, inclusive end) from which to read data");
                    endColField.setText("2147483646");
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
                ((GridBagLayout)buttonBar.getLayout()).columnWidths = new int[] {0, 85, 80};
                ((GridBagLayout)buttonBar.getLayout()).columnWeights = new double[] {1.0, 0.0, 0.0};

                //---- saveButton ----
                saveButton.setText("Save Profile");
                saveButton.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        saveButtonActionPerformed(e);
                    }
                });
                buttonBar.add(saveButton, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 0, 5), 0, 0));

                //---- okButton ----
                okButton.setText("Load Data");
                okButton.setVisible(false);
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
        pack();
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JPanel panel16;
    private JLabel label16;
    private JComboBox DBSubprotocol;
    private JPanel panel4;
    private JLabel label4;
    private JTextField DBName;
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
    private JSeparator separator1;
    private JPanel panel15;
    private JLabel label15;
    private JTextField newTableNameField;
    private JButton addNewTableButton;
    private JSeparator separator2;
    private JPanel panel7;
    private JLabel label7;
    private JComboBox tableName;
    private JPanel panel8;
    private JLabel label8;
    private JComboBox dataType;
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


    private JComponent[] tableFields;

}
