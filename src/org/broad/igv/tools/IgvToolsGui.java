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

package org.broad.igv.tools;

import org.apache.log4j.Logger;
import org.broad.igv.sam.SpliceJunctionFinderTrack;
import org.broad.igv.track.WindowFunction;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import javax.swing.border.TitledBorder;
import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;

public class IgvToolsGui extends JDialog {

    private static Logger log = Logger.getLogger(IgvToolsGui.class);

    static JFileChooser fileDialog;

    public enum Tool {
        COUNT("Count"),
        SORT("Sort"),
        INDEX("Index"),
        TILE("toTDF");

        private String displayName;

        Tool(String displayName) {
            this.displayName = displayName;
        }

        @Override
        public String toString() {
            return this.displayName;
        }
    }

    String[] zoomLevels = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"};

    PrintStream systemOutStream;
    PrintStream systemErrStream;

    static File lastDirectory;

    IgvTools igvTools;

    public IgvToolsGui() {

        this.igvTools = new IgvTools();
        initComponents();
        initUI();
        updateUI();

        closeButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                close();
                setVisible(false);
            }
        });
        inputButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                try {
                    File chosenFile = chooseFile();
                    inputField.setText(chosenFile.getAbsolutePath());
                    updateUI();
                } catch (NullPointerException npe) {
                }
            }
        });
        outputButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                try {
                    File chosenFile = chooseFile();
                    outputField.setText(chosenFile.getAbsolutePath());
                    updateUI();
                } catch (NullPointerException npe) {
                }
            }
        });
        genomeButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                try {
                    File chosenFile = chooseFile();
                    genomeField.setText(chosenFile.getAbsolutePath());
                    updateUI();
                } catch (NullPointerException npe) {
                }
            }
        });
        probeButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                try {
                    File chosenFile = chooseFile();
                    probeField.setText(chosenFile.getAbsolutePath());
                    updateUI();
                } catch (NullPointerException npe) {
                }
            }
        });
        toolCombo.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                updateUI();
            }
        });
        runButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                run();
            }
        });

        tempButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                try {
                    File chosenFile = chooseFile();
                    tmpDirectoryField.setText(chosenFile.getAbsolutePath());
                    updateUI();
                } catch (NullPointerException npe) {
                }
            }
        });
    }

    private void initUI() {
        setContentPane(mainPanel);
        setModal(true);

        addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent windowEvent) {
                close();
            }
        });

        for (String item : zoomLevels) {
            zoomCombo.addItem(item);
        }
        for (Tool tool : Tool.values()) {
            toolCombo.addItem(tool);
        }

        zoomCombo.setSelectedIndex(IgvTools.MAX_ZOOM);
        windowSizeField.setText(String.valueOf(IgvTools.WINDOW_SIZE));
        maxRecordsField.setText(String.valueOf(IgvTools.MAX_RECORDS_IN_RAM));
        redirectSystemStreams();
    }

    private void close() {
        System.setErr(systemErrStream);
        System.setOut(systemOutStream);
        dispose();
    }

    private void updateTextArea(final String text) {
        SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                outputText.append(text);
            }
        });
    }

    private void redirectSystemStreams() {
        OutputStream out = new OutputStream() {
            @Override
            public void write(int b) throws IOException {
                updateTextArea(String.valueOf((char) b));
            }

            @Override
            public void write(byte[] b, int off, int len) throws IOException {
                updateTextArea(new String(b, off, len));
            }

            @Override
            public void write(byte[] b) throws IOException {
                write(b, 0, b.length);
            }
        };

        systemOutStream = System.out;
        systemErrStream = System.err;

        System.setOut(new PrintStream(out, true));
        System.setErr(new PrintStream(out, true));
    }


    private void updateUI() {
        Tool tool = (Tool) toolCombo.getSelectedItem();

        if (tool.equals(Tool.COUNT)) {
            inputField.setEnabled(true);
            inputButton.setEnabled(true);
            outputField.setEnabled(true);
            outputButton.setEnabled(true);
            outputLabel.setEnabled(true);

            if (!this.genomeSelectionDisabled) {
                genomeField.setEnabled(true);
                genomeButton.setEnabled(true);
                genomeLabel.setEnabled(true);
            }

            zoomCombo.setEnabled(true);
            zoomLabel.setEnabled(true);
            maxRecordsField.setEnabled(false);
            maxRecordsLabel.setEnabled(false);
            probeField.setEnabled(false);
            probeButton.setEnabled(false);
            probeLabel.setEnabled(false);
            tmpDirectoryLabel.setEnabled(false);
            tmpDirectoryField.setEnabled(false);
            tempButton.setEnabled(false);
            windowSizeLabel.setEnabled(true);
            windowSizeField.setEnabled(true);

            pairsCB.setEnabled(true);
            extFactorLabel.setEnabled(true);
            extFactorField.setEnabled(true);

            enableWindowFunctions();

        } else if (tool.equals(Tool.SORT)) {
            inputField.setEnabled(true);
            inputButton.setEnabled(true);
            outputField.setEnabled(true);
            outputButton.setEnabled(true);
            outputLabel.setEnabled(true);
            genomeField.setEnabled(false);
            genomeButton.setEnabled(false);
            genomeLabel.setEnabled(false);
            zoomCombo.setEnabled(false);
            zoomLabel.setEnabled(false);
            maxRecordsField.setEnabled(true);
            maxRecordsLabel.setEnabled(true);
            probeField.setEnabled(false);
            probeButton.setEnabled(false);
            probeLabel.setEnabled(false);
            tmpDirectoryLabel.setEnabled(true);
            tmpDirectoryField.setEnabled(true);
            tempButton.setEnabled(true);
            windowSizeLabel.setEnabled(false);
            windowSizeField.setEnabled(false);
            pairsCB.setEnabled(false);
            extFactorLabel.setEnabled(false);
            extFactorField.setEnabled(false);
            disableWindowFunctions();

        } else if (tool.equals(Tool.INDEX)) {
            inputField.setEnabled(true);
            inputButton.setEnabled(true);
            outputField.setEnabled(false);
            outputButton.setEnabled(false);
            outputLabel.setEnabled(false);
            genomeField.setEnabled(false);
            genomeButton.setEnabled(false);
            genomeLabel.setEnabled(false);
            zoomCombo.setEnabled(false);
            zoomLabel.setEnabled(false);
            maxRecordsField.setEnabled(false);
            maxRecordsLabel.setEnabled(false);
            probeField.setEnabled(false);
            probeButton.setEnabled(false);
            probeLabel.setEnabled(false);
            tmpDirectoryLabel.setEnabled(false);
            tmpDirectoryField.setEnabled(false);
            tempButton.setEnabled(false);
            windowSizeLabel.setEnabled(false);
            windowSizeField.setEnabled(false);
            pairsCB.setEnabled(false);
            extFactorLabel.setEnabled(false);
            extFactorField.setEnabled(false);
            disableWindowFunctions();

        } else if (tool.equals(Tool.TILE)) {
            inputField.setEnabled(true);
            inputButton.setEnabled(true);
            outputField.setEnabled(true);
            outputButton.setEnabled(true);
            outputLabel.setEnabled(true);

            if (!this.genomeSelectionDisabled) {
                genomeField.setEnabled(true);
                genomeButton.setEnabled(true);
                genomeLabel.setEnabled(true);
            }

            zoomCombo.setEnabled(true);
            zoomLabel.setEnabled(true);
            maxRecordsField.setEnabled(true);
            maxRecordsLabel.setEnabled(true);
            probeField.setEnabled(true);
            probeButton.setEnabled(true);
            probeLabel.setEnabled(true);
            tmpDirectoryLabel.setEnabled(true);
            tmpDirectoryField.setEnabled(true);
            tempButton.setEnabled(false);
            windowSizeLabel.setEnabled(false);
            windowSizeField.setEnabled(false);
            pairsCB.setEnabled(false);
            extFactorLabel.setEnabled(false);
            extFactorField.setEnabled(false);
            enableWindowFunctions();
        }
    }

    private void enableWindowFunctions() {
        Component[] com = windowFunctionPanel.getComponents();
        windowFunctionLabel.setEnabled(true);
        for (int i = 0; i < com.length; i++) {
            com[i].setEnabled(true);
        }
    }

    private void disableWindowFunctions() {
        Component[] com = windowFunctionPanel.getComponents();
        windowFunctionLabel.setEnabled(false);
        for (int i = 0; i < com.length; i++) {
            com[i].setEnabled(false);
        }
    }


    private boolean validateFields() {
        Tool tool = (Tool) toolCombo.getSelectedItem();

        // All commands require in input file
        if (inputField.getText().trim().length() == 0) {
            showMessage("Input file is required");
            return false;
        }

        if (tool.equals(Tool.INDEX)) {
            return true;
        }

        // All commands except "index" require an output file
        if (outputField.getText().trim().length() == 0) {
            showMessage("Output file is required");
            return false;
        }

        // See if file exists
        File outputFile = new File(outputField.getText().trim());
        if (outputFile.exists()) {
            int opt = JOptionPane.showConfirmDialog(this, "Output file: " + outputField.getText() + " exists.  Overwite?");
            if (opt != JOptionPane.YES_OPTION) {
                return false;
            }
        }

        //Check that parent directory exists
        if (!outputFile.getParentFile().exists()) {
            showMessage("Directory " + outputFile.getParent() + " does not exist. Output must use an existing directory");
            return false;
        }

        if (tool.equals(Tool.SORT)) {
            return true;
        }

        if (genomeField.getText().trim().length() == 0) {
            showMessage("Genome is required");
            return false;
        }

        return true;
    }

    private void run() {
        Tool tool = (Tool) toolCombo.getSelectedItem();
        if (!validateFields()) {
            return;
        }
        try {
            switch (tool) {
                case COUNT:
                    doCount();
                    break;
                case SORT:
                    doSort();
                    break;
                case INDEX:
                    doIndex();
                    break;
                case TILE:
                    doTile();
                    break;
            }
        } catch (PreprocessingException e) {
            showMessage("Error performing " + tool + ": " + e.getMessage());
        }
    }

    private void showMessage(String tool) {
        JOptionPane.showMessageDialog(this, tool);
    }

    private void doSort() {
        // Sort is not interruptible, and there is no provision for progress.  So use a thread & wait cursor
        runButton.setEnabled(false);
        SwingWorker swingWorker = new IgvToolsSwingWorker() {

            @Override
            protected Object doInBackground() {
                try {
                    setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
                    String maxRecordText = maxRecordsField.getText();
                    int maxRecords = (maxRecordText != null && maxRecordText.length() > 0) ?
                            Integer.parseInt(maxRecordText) : IgvTools.MAX_RECORDS_IN_RAM;
                    igvTools.doSort(inputField.getText(), outputField.getText(), tmpDirectoryField.getText(), maxRecords);
                } catch (Exception e) {
                    showMessage("Error: " + e.getMessage());
                }

                return null;
            }
        };

        swingWorker.execute();
    }

    private void doCount() {

        SwingWorker swingWorker = new IgvToolsSwingWorker() {

            @Override
            protected Object doInBackground() {
                try {
                    setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
                    String ifile = inputField.getText();
                    String ofile = outputField.getText();
                    String genomeId = genomeField.getText();
                    int maxZoomValue = Integer.parseInt(zoomCombo.getSelectedItem().toString());
                    Collection<WindowFunction> wfs = getWindowFunctions();

                    String windowSizeText = windowSizeField.getText();
                    int windowSize = (windowSizeText != null && windowSizeText.length() > 0) ?
                            Integer.parseInt(windowSizeText) : IgvTools.WINDOW_SIZE;

                    String extFactorText = extFactorField.getText().trim();
                    int extFactor = extFactorText.length() == 0 ? 0 : Integer.parseInt(extFactorText);
                    int strandOption = pairsCB.isSelected() ? CoverageCounter.PAIRED_COVERAGE : 0;

                    runButton.setEnabled(false);
                    int preExtFactor = 0;
                    int postExtFactor = 0;
                    igvTools.doCount(ifile, ofile, genomeId, maxZoomValue, wfs, windowSize, extFactor,
                            preExtFactor, postExtFactor, null, null, 0, strandOption);
                } catch (Exception e) {
                    log.error(e.getMessage(), e);
                    showMessage("Error: " + e.getMessage());
                }

                return null;
            }
        };

        swingWorker.execute();
    }

    private void doTile() {

        SwingWorker swingWorker = new IgvToolsSwingWorker() {

            @Override
            protected Object doInBackground() {
                try {
                    setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
                    String ifile = inputField.getText();
                    String ofile = outputField.getText();
                    String genomeId = genomeField.getText();
                    int maxZoomValue = Integer.parseInt(zoomCombo.getSelectedItem().toString());
                    Collection<WindowFunction> wfs = getWindowFunctions();
                    String probeFile = probeField.getText();

                    String typeString = Preprocessor.getExtension(ifile);

                    runButton.setEnabled(false);
                    igvTools.toTDF(typeString, ifile, ofile, probeFile, genomeId, maxZoomValue, wfs, null, IgvTools.MAX_RECORDS_IN_RAM);
                } catch (Exception e) {
                    showMessage("Error: " + e.getMessage());
                }

                return null;
            }
        };
        swingWorker.execute();
    }

    //    public static void doIndex(String ifile, int indexType, int binSize) throws IOException {

    private void doIndex() {

        SwingWorker swingWorker = new IgvToolsSwingWorker() {

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
                    log.error("Error indexing file", e);
                    showMessage("Error: " + e.getMessage());
                }

                return null;
            }
        };

        swingWorker.execute();
    }

    private Collection<WindowFunction> getWindowFunctions() {
        ArrayList<WindowFunction> wfs = new ArrayList();

        if (minCheckBox.isSelected()) {
            wfs.add(WindowFunction.min);
        }
        if (maxCheckBox.isSelected()) {
            wfs.add(WindowFunction.max);
        }
        if (meanCheckBox.isSelected()) {
            wfs.add(WindowFunction.mean);
        }
        if (a98CheckBox.isSelected()) {
            wfs.add(WindowFunction.percentile98);
        }
        if (a90CheckBox.isSelected()) {
            wfs.add(WindowFunction.percentile90);
        }
        if (a10CheckBox.isSelected()) {
            wfs.add(WindowFunction.percentile10);
        }
        if (a2CheckBox.isSelected()) {
            wfs.add(WindowFunction.percentile2);
        }

        if (wfs.isEmpty()) {
            wfs.add(WindowFunction.mean);
        }

        return wfs;
    }

    /**
     * Select an input file, or directory in some cases
     *
     * @return
     */
    private File chooseFile() {

        //TODO Override so you can specify the file type with a string array ex: {".wig", ".tdf"}

        // if (fileDialog == null) {
        fileDialog = new JFileChooser();
        // }

        if (lastDirectory != null) {
            fileDialog.setCurrentDirectory(lastDirectory);
        }

        fileDialog.setMultiSelectionEnabled(false);

        fileDialog.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);

        int returnVal = fileDialog.showDialog(this, "Select File");
        if (returnVal == JFileChooser.CANCEL_OPTION) {
            return null;
        } else {
            File selected = fileDialog.getSelectedFile();
            lastDirectory = selected.isDirectory() ? selected : selected.getParentFile();
            return selected;
        }
    }

    public static void main(String[] args) {
        launch(true, null);
    }

    private boolean genomeSelectionDisabled = false;

    public static void launch(boolean modal, String genomeId) {
        IgvToolsGui mainWindow = new IgvToolsGui();
        mainWindow.pack();
        mainWindow.setModal(modal);
        mainWindow.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        mainWindow.setResizable(false);

        if (genomeId != null) {
            mainWindow.genomeField.setText(genomeId);
            mainWindow.genomeField.setEnabled(false);
            mainWindow.genomeField.setToolTipText("<html>To change the genome id close this window and <br>use the pulldown on the IGV batch screen.");
            mainWindow.genomeButton.setEnabled(false);
            mainWindow.genomeSelectionDisabled = true;
        }

        mainWindow.setVisible(true);
    }


    private void inputButtonActionPerformed(ActionEvent e) {
        setDefaultOutputText();
    }

    private void inputFieldActionPerformed(ActionEvent e) {
        setDefaultOutputText();
    }

    private void inputFieldFocusLost(FocusEvent e) {
        setDefaultOutputText();
    }

    /**
     * Return what the defaultOutputText should be. {@code null} implies
     * there should be no change
     *
     * @param inputFieldText
     * @param tool
     * @return
     */
    static String getDefaultOutputText(String inputFieldText, Tool tool) {
        String defaultOutputText = null;
        if (inputFieldText.length() > 0) {
            switch (tool) {
                case COUNT:
                case TILE:
                    defaultOutputText = inputFieldText + ".tdf";
                    break;
                case SORT:
                    int ext = inputFieldText.lastIndexOf(".");
                    if (ext > 0) {
                        defaultOutputText = inputFieldText.substring(0, ext) + ".sorted" + inputFieldText.substring(ext);
                    }
                    break;
                case INDEX:
                    defaultOutputText = "";
                    break;
            }
        }
        return defaultOutputText;
    }

    private void setDefaultOutputText() {
        String defaultOutputText = getDefaultOutputText(inputField.getText(), (Tool) toolCombo.getSelectedItem());
        if (defaultOutputText != null) outputField.setText(defaultOutputText);
    }

    private void toolComboActionPerformed(ActionEvent e) {
        setDefaultOutputText();
    }


    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        mainPanel = new JPanel();
        requiredPanel = new JPanel();
        toolCombo = new JComboBox();
        JLabel label1 = new JLabel();
        JLabel label2 = new JLabel();
        outputLabel = new JLabel();
        outputButton = new JButton();
        inputField = new JTextField();
        inputButton = new JButton();
        outputField = new JTextField();
        genomeLabel = new JLabel();
        genomeField = new JTextField();
        genomeButton = new JButton();
        tilePanel = new JPanel();
        zoomLabel = new JLabel();
        windowFunctionLabel = new JLabel();
        windowFunctionPanel = new JPanel();
        minCheckBox = new JCheckBox();
        maxCheckBox = new JCheckBox();
        meanCheckBox = new JCheckBox();
        medianCheckBox = new JCheckBox();
        a2CheckBox = new JCheckBox();
        a10CheckBox = new JCheckBox();
        a90CheckBox = new JCheckBox();
        a98CheckBox = new JCheckBox();
        probeLabel = new JLabel();
        probeField = new JTextField();
        probeButton = new JButton();
        zoomCombo = new JComboBox();
        windowSizeLabel = new JLabel();
        windowSizeField = new JTextField();
        extFactorLabel = new JLabel();
        extFactorField = new JTextField();
        pairsCB = new JCheckBox();
        sortPanel = new JPanel();
        tmpDirectoryLabel = new JLabel();
        maxRecordsLabel = new JLabel();
        tempButton = new JButton();
        tmpDirectoryField = new JTextField();
        maxRecordsField = new JTextField();
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

                //---- toolCombo ----
                toolCombo.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        toolComboActionPerformed(e);
                    }
                });
                requiredPanel.add(toolCombo, new GridBagConstraints(2, 1, 1, 1, 1.0, 1.0,
                    GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL,
                    new Insets(0, 0, 0, 0), 0, 0));

                //---- label1 ----
                label1.setText("Command");
                requiredPanel.add(label1, new GridBagConstraints(1, 1, 1, 1, 0.0, 1.0,
                    GridBagConstraints.WEST, GridBagConstraints.NONE,
                    new Insets(0, 0, 0, 0), 0, 0));

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
                    @Override
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
                    @Override
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

                //---- genomeLabel ----
                genomeLabel.setToolTipText("Either a genome ID (e.g. hg18) or the full path to a .genome file.");
                genomeLabel.setText("Genome");
                requiredPanel.add(genomeLabel, new GridBagConstraints(1, 4, 1, 1, 0.0, 1.0,
                    GridBagConstraints.WEST, GridBagConstraints.NONE,
                    new Insets(0, 0, 0, 0), 0, 0));
                requiredPanel.add(genomeField, new GridBagConstraints(2, 4, 1, 1, 1.0, 1.0,
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

            //======== tilePanel ========
            {
                tilePanel.setEnabled(true);
                tilePanel.setFont(tilePanel.getFont().deriveFont(Font.ITALIC, 10f));
                tilePanel.setBorder(new TitledBorder(null, "TDF and Count options", TitledBorder.LEADING, TitledBorder.TOP));
                tilePanel.setLayout(new GridBagLayout());

                //---- zoomLabel ----
                zoomLabel.setToolTipText("<html>Specifies the maximum zoom level to precompute. The default value is 7.<br>To reduce file size at the expense of Iperformance this value can be reduced.");
                zoomLabel.setText("Zoom Levels");
                tilePanel.add(zoomLabel, new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0,
                    GridBagConstraints.WEST, GridBagConstraints.NONE,
                    new Insets(0, 0, 0, 0), 0, 0));

                //---- windowFunctionLabel ----
                windowFunctionLabel.setToolTipText("Window functions to use for summarizing data. ");
                windowFunctionLabel.setText("Window Functions");
                tilePanel.add(windowFunctionLabel, new GridBagConstraints(1, 2, 1, 1, 0.0, 0.0,
                    GridBagConstraints.NORTHWEST, GridBagConstraints.NONE,
                    new Insets(0, 0, 0, 0), 0, 0));

                //======== windowFunctionPanel ========
                {
                    windowFunctionPanel.setLayout(new GridLayout(2, 0));

                    //---- minCheckBox ----
                    minCheckBox.setText("Min");
                    windowFunctionPanel.add(minCheckBox);

                    //---- maxCheckBox ----
                    maxCheckBox.setText("Max");
                    windowFunctionPanel.add(maxCheckBox);

                    //---- meanCheckBox ----
                    meanCheckBox.setSelected(true);
                    meanCheckBox.setText("Mean");
                    windowFunctionPanel.add(meanCheckBox);

                    //---- medianCheckBox ----
                    medianCheckBox.setText("Median");
                    windowFunctionPanel.add(medianCheckBox);

                    //---- a2CheckBox ----
                    a2CheckBox.setText("2%");
                    windowFunctionPanel.add(a2CheckBox);

                    //---- a10CheckBox ----
                    a10CheckBox.setText("10%");
                    windowFunctionPanel.add(a10CheckBox);

                    //---- a90CheckBox ----
                    a90CheckBox.setText("90%");
                    windowFunctionPanel.add(a90CheckBox);

                    //---- a98CheckBox ----
                    a98CheckBox.setText("98%");
                    windowFunctionPanel.add(a98CheckBox);
                }
                tilePanel.add(windowFunctionPanel, new GridBagConstraints(2, 2, 1, 1, 1.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 0, 0), 0, 0));

                //---- probeLabel ----
                probeLabel.setFont(probeLabel.getFont());
                probeLabel.setToolTipText("<html>Specifies a \"bed\" file to be used to map probe identifiers to locations.  This option is useful <br>when preprocessing gct files.  The bed file should contain 4 columns: chr start end name\n<br>where name is the probe name in the gct file.");
                probeLabel.setText("Probe to Loci Mapping");
                tilePanel.add(probeLabel, new GridBagConstraints(1, 3, 1, 1, 0.0, 1.0,
                    GridBagConstraints.WEST, GridBagConstraints.NONE,
                    new Insets(0, 0, 0, 0), 0, 0));
                tilePanel.add(probeField, new GridBagConstraints(2, 3, 1, 1, 1.0, 1.0,
                    GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL,
                    new Insets(0, 0, 0, 0), 0, 0));

                //---- probeButton ----
                probeButton.setText("Browse");
                tilePanel.add(probeButton, new GridBagConstraints(3, 3, 1, 1, 0.0, 1.0,
                    GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL,
                    new Insets(0, 0, 0, 0), 0, 0));

                //---- zoomCombo ----
                zoomCombo.setEditable(false);
                zoomCombo.setModel(new DefaultComboBoxModel(new String[] {

                }));
                tilePanel.add(zoomCombo, new GridBagConstraints(2, 1, 1, 1, 1.0, 0.0,
                    GridBagConstraints.WEST, GridBagConstraints.NONE,
                    new Insets(0, 0, 0, 0), 0, 0));

                //---- windowSizeLabel ----
                windowSizeLabel.setToolTipText("The window size over which coverage computed when using the count command.  Defaults to 25 bp.");
                windowSizeLabel.setText("Window Size");
                tilePanel.add(windowSizeLabel, new GridBagConstraints(1, 4, 1, 1, 0.0, 0.0,
                    GridBagConstraints.WEST, GridBagConstraints.NONE,
                    new Insets(0, 0, 0, 0), 0, 0));
                tilePanel.add(windowSizeField, new GridBagConstraints(2, 4, 1, 1, 1.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL,
                    new Insets(0, 0, 0, 0), 0, 0));

                //---- extFactorLabel ----
                extFactorLabel.setText("Extension Factor");
                extFactorLabel.setToolTipText("Extend read by this amount.  Set to average fragment size of library");
                tilePanel.add(extFactorLabel, new GridBagConstraints(1, 5, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 0, 0), 0, 0));
                tilePanel.add(extFactorField, new GridBagConstraints(2, 5, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 0, 0), 0, 0));

                //---- pairsCB ----
                pairsCB.setText("Count as Pairs");
                pairsCB.setToolTipText("Count area between proper pairs as covered.");
                tilePanel.add(pairsCB, new GridBagConstraints(1, 6, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 0, 0), 0, 0));
            }
            mainPanel.add(tilePanel, new GridBagConstraints(1, 2, 1, 1, 1.0, 1.0,
                GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                new Insets(0, 0, 10, 0), 0, 0));

            //======== sortPanel ========
            {
                sortPanel.setBorder(new TitledBorder(null, "Sort Options", TitledBorder.LEADING, TitledBorder.TOP));
                sortPanel.setLayout(new GridBagLayout());

                //---- tmpDirectoryLabel ----
                tmpDirectoryLabel.setToolTipText("<html>Specify a temporary working directory.  For large input files this directory will be used to <br>store intermediate results of the sort. The default is the users temp directory.");
                tmpDirectoryLabel.setText("Temp Directory");
                sortPanel.add(tmpDirectoryLabel, new GridBagConstraints(1, 1, 1, 1, 0.0, 1.0,
                    GridBagConstraints.WEST, GridBagConstraints.NONE,
                    new Insets(0, 0, 0, 0), 0, 0));

                //---- maxRecordsLabel ----
                maxRecordsLabel.setToolTipText("<html>The maximum number of records to keep in memory during the sort.  The default value is <br>500000.  Increase this number if you receive \"too many open files\" errors.   Decrease it if you <br>experience \"out of memory\" errors.");
                maxRecordsLabel.setText("Max Records");
                sortPanel.add(maxRecordsLabel, new GridBagConstraints(1, 2, 1, 1, 0.0, 1.0,
                    GridBagConstraints.WEST, GridBagConstraints.NONE,
                    new Insets(0, 0, 0, 0), 0, 0));

                //---- tempButton ----
                tempButton.setText("Browse");
                sortPanel.add(tempButton, new GridBagConstraints(4, 1, 1, 1, 0.0, 1.0,
                    GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL,
                    new Insets(0, 0, 0, 0), 0, 0));
                sortPanel.add(tmpDirectoryField, new GridBagConstraints(3, 1, 1, 1, 1.0, 1.0,
                    GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL,
                    new Insets(0, 0, 0, 0), 0, 0));
                sortPanel.add(maxRecordsField, new GridBagConstraints(3, 2, 1, 1, 1.0, 1.0,
                    GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL,
                    new Insets(0, 0, 0, 0), 0, 0));
            }
            mainPanel.add(sortPanel, new GridBagConstraints(1, 3, 1, 1, 1.0, 1.0,
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
                    for(int i = 0; i < OutputPanel.getComponentCount(); i++) {
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
    private JComboBox toolCombo;
    private JLabel outputLabel;
    private JButton outputButton;
    private JTextField inputField;
    private JButton inputButton;
    private JTextField outputField;
    private JLabel genomeLabel;
    private JTextField genomeField;
    private JButton genomeButton;
    private JPanel tilePanel;
    private JLabel zoomLabel;
    private JLabel windowFunctionLabel;
    private JPanel windowFunctionPanel;
    private JCheckBox minCheckBox;
    private JCheckBox maxCheckBox;
    private JCheckBox meanCheckBox;
    private JCheckBox medianCheckBox;
    private JCheckBox a2CheckBox;
    private JCheckBox a10CheckBox;
    private JCheckBox a90CheckBox;
    private JCheckBox a98CheckBox;
    private JLabel probeLabel;
    private JTextField probeField;
    private JButton probeButton;
    private JComboBox zoomCombo;
    private JLabel windowSizeLabel;
    private JTextField windowSizeField;
    private JLabel extFactorLabel;
    private JTextField extFactorField;
    private JCheckBox pairsCB;
    private JPanel sortPanel;
    private JLabel tmpDirectoryLabel;
    private JLabel maxRecordsLabel;
    private JButton tempButton;
    private JTextField tmpDirectoryField;
    private JTextField maxRecordsField;
    private JPanel buttonPanel;
    private JButton runButton;
    private JButton closeButton;
    private JPanel OutputPanel;
    private JScrollPane outputScroll;
    private JTextArea outputText;
    private JProgressBar progressBar;
    // JFormDesigner - End of variables declaration  //GEN-END:variables

    private abstract class IgvToolsSwingWorker extends SwingWorker{

        @Override
        protected void done() {
            runButton.setEnabled(true);
            setCursor(Cursor.getDefaultCursor());
            updateTextArea("Done");
        }
    }

}
