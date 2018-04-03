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
 * Created by JFormDesigner on Tue Nov 09 10:32:35 EST 2010
 */

package org.broad.igv.tools.ui;

import org.broad.igv.tools.IgvTools;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import javax.swing.border.TitledBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;

/**
 * @author Stan Diamond
 */
public class IndexGui extends JDialog {

    IgvTools igvTools;

    public IndexGui() {

        this.igvTools = new IgvTools();
        initComponents();

        setTitle("Create Index");

        // Hide output file field for now
        this.outputLabel.setVisible(true);
        this.outputField.setVisible(true);
        this.outputButton.setVisible(true);
    }

    private void inputFieldFocusLost(FocusEvent e) {
        // TODO add your code here
    }

    private void inputFieldActionPerformed(ActionEvent e) {
        // TODO add your code here
    }

    private void inputButtonActionPerformed(ActionEvent e) {
        // TODO add your code here
    }

    //    public static void doIndex(String ifile, int indexType, int binSize) throws IOException {

    private void doIndex() {

        SwingWorker swingWorker = new SwingWorker() {

            @Override
            protected Object doInBackground() {
                try {
                    setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
                    String ifile = inputField.getText();
                    int indexType = IgvTools.LINEAR_INDEX;
                    int binSize = IgvTools.LINEAR_BIN_SIZE;

                    runButton.setEnabled(false);
                    igvTools.doIndex(ifile, null, indexType, binSize);
                } catch (Exception e) {
                    showMessage("Error: " + e.getMessage());
                }

                return null;
            }

            @Override
            protected void done() {
                runButton.setEnabled(true);
                setCursor(Cursor.getDefaultCursor());
            }
        };

        swingWorker.execute();
    }


    private void showMessage(String tool) {
        JOptionPane.showMessageDialog(this, tool);
    }


    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        mainPanel = new JPanel();
        requiredPanel = new JPanel();
        JLabel label2 = new JLabel();
        outputLabel = new JLabel();
        outputButton = new JButton();
        inputField = new JTextField();
        inputButton = new JButton();
        outputField = new JTextField();
        genomeButton = new JButton();
        buttonPanel = new JPanel();
        runButton = new JButton();
        closeButton = new JButton();
        OutputPanel = new JPanel();
        outputScroll = new JScrollPane();
        outputText = new JTextArea();
        JSeparator separator1 = new JSeparator();
        progressBar = new JProgressBar();

        //======== mainPanel ========
        {
            mainPanel.setBorder(new TitledBorder(new EtchedBorder(), ""));
            mainPanel.setLayout(new GridBagLayout());

            //======== requiredPanel ========
            {
                requiredPanel.setLayout(new GridBagLayout());

                //---- label2 ----
                label2.setText("Input File");
                requiredPanel.add(label2, new GridBagConstraints(1, 2, 1, 1, 0.0, 1.0,
                        GridBagConstraints.WEST, GridBagConstraints.NONE,
                        new Insets(0, 0, 0, 0), 0, 0));

                //---- outputLabel ----
                outputLabel.setText("Output File");
                requiredPanel.add(outputLabel, new GridBagConstraints(1, 3, 1, 1, 0.0, 1.0,
                        GridBagConstraints.WEST, GridBagConstraints.NONE,
                        new Insets(0, 0, 0, 0), 0, 0));

                //---- outputButton ----
                outputButton.setText("Browse");
                requiredPanel.add(outputButton, new GridBagConstraints(3, 3, 1, 1, 0.0, 1.0,
                        GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL,
                        new Insets(0, 0, 0, 0), 0, 0));

                //---- inputField ----
                inputField.addFocusListener(new FocusAdapter() {
                    @Override
                    public void focusLost(FocusEvent e) {
                        inputFieldFocusLost(e);
                    }
                });
                inputField.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        inputFieldActionPerformed(e);
                    }
                });
                requiredPanel.add(inputField, new GridBagConstraints(2, 2, 1, 1, 1.0, 1.0,
                        GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL,
                        new Insets(0, 0, 0, 0), 0, 0));

                //---- inputButton ----
                inputButton.setText("Browse");
                inputButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        inputButtonActionPerformed(e);
                    }
                });
                requiredPanel.add(inputButton, new GridBagConstraints(3, 2, 1, 1, 0.0, 1.0,
                        GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL,
                        new Insets(0, 0, 0, 0), 0, 0));
                requiredPanel.add(outputField, new GridBagConstraints(2, 3, 1, 1, 1.0, 1.0,
                        GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL,
                        new Insets(0, 0, 0, 0), 0, 0));

                //---- genomeButton ----
                genomeButton.setText("Browse");
                requiredPanel.add(genomeButton, new GridBagConstraints(3, 4, 1, 1, 0.0, 1.0,
                        GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL,
                        new Insets(0, 0, 0, 0), 0, 0));
            }
            mainPanel.add(requiredPanel, new GridBagConstraints(1, 1, 1, 1, 1.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 10, 0), 0, 0));

            //======== buttonPanel ========
            {
                buttonPanel.setLayout(new GridBagLayout());

                //---- runButton ----
                runButton.setText("Run");
                buttonPanel.add(runButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 1.0,
                        GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL,
                        new Insets(0, 0, 0, 0), 0, 0));

                //---- closeButton ----
                closeButton.setText("Close");
                buttonPanel.add(closeButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 1.0,
                        GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL,
                        new Insets(0, 0, 0, 0), 0, 0));
            }
            mainPanel.add(buttonPanel, new GridBagConstraints(1, 4, 1, 1, 1.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 10, 0), 0, 0));

            //======== OutputPanel ========
            {
                OutputPanel.setBorder(new TitledBorder(BorderFactory.createEmptyBorder(), "Messages", TitledBorder.LEADING, TitledBorder.TOP));
                OutputPanel.setLayout(null);

                //======== outputScroll ========
                {

                    //---- outputText ----
                    outputText.setEditable(false);
                    outputText.setText("");
                    outputText.setRows(10);
                    outputScroll.setViewportView(outputText);
                }
                OutputPanel.add(outputScroll);
                outputScroll.setBounds(4, 20, 881, outputScroll.getPreferredSize().height);

                { // compute preferred size
                    Dimension preferredSize = new Dimension();
                    for (int i = 0; i < OutputPanel.getComponentCount(); i++) {
                        Rectangle bounds = OutputPanel.getComponent(i).getBounds();
                        preferredSize.width = Math.max(bounds.x + bounds.width, preferredSize.width);
                        preferredSize.height = Math.max(bounds.y + bounds.height, preferredSize.height);
                    }
                    Insets insets = OutputPanel.getInsets();
                    preferredSize.width += insets.right;
                    preferredSize.height += insets.bottom;
                    OutputPanel.setMinimumSize(preferredSize);
                    OutputPanel.setPreferredSize(preferredSize);
                }
            }
            mainPanel.add(OutputPanel, new GridBagConstraints(1, 6, 1, 1, 1.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 10, 0), 0, 0));
            mainPanel.add(separator1, new GridBagConstraints(1, 5, 1, 1, 1.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 10, 0), 0, 0));
            mainPanel.add(progressBar, new GridBagConstraints(1, 7, 1, 1, 1.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL,
                    new Insets(0, 0, 0, 0), 0, 0));
        }
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel mainPanel;
    private JPanel requiredPanel;
    private JLabel outputLabel;
    private JButton outputButton;
    private JTextField inputField;
    private JButton inputButton;
    private JTextField outputField;
    private JButton genomeButton;
    private JPanel buttonPanel;
    private JButton runButton;
    private JButton closeButton;
    private JPanel OutputPanel;
    private JScrollPane outputScroll;
    private JTextArea outputText;
    private JProgressBar progressBar;
    // JFormDesigner - End of variables declaration  //GEN-END:variables


    public static void launch(boolean modal) {
        IndexGui mainWindow = new IndexGui();
        mainWindow.pack();
        mainWindow.setModal(modal);
        mainWindow.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        mainWindow.setResizable(false);


        mainWindow.setVisible(true);
    }
}
