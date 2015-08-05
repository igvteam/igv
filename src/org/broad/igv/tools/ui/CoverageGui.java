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
 * Created by JFormDesigner on Mon Nov 08 21:24:47 EST 2010
 */

package org.broad.igv.tools.ui;

import org.broad.igv.PreferenceManager;
import org.broad.igv.tools.CoverageCounter;
import org.broad.igv.tools.IgvTools;
import org.broad.igv.tools.Preprocessor;
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

/**
 * @author Stan Diamond
 */
public class CoverageGui extends JDialog {

    public enum Mode {
        COVERAGE, TILE
    }

    Mode mode;
    String[] zoomLevels = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"};
    PrintStream systemOutStream;
    PrintStream systemErrStream;
    IgvTools igvTools;
    private JFileChooser fileDialog;

    public CoverageGui(Mode mode) {
        super();
        this.mode = mode;
        initComponents();
        initUI();
    }


    private void inputButtonActionPerformed(ActionEvent e) {
        File chosenFile = chooseFile();
        inputField.setText(chosenFile.getAbsolutePath());
        setDefaultOutputText();
    }

    private void inputFieldActionPerformed(ActionEvent e) {
        setDefaultOutputText();
    }

    private void inputFieldFocusLost(FocusEvent e) {
        setDefaultOutputText();
    }

    private void setDefaultOutputText() {
        if (inputField.getText().length() > 0) {
            outputField.setText(inputField.getText() + ".tdf");
        }
    }

    private void runButtonActionPerformed(ActionEvent e) {
        if (mode == Mode.COVERAGE) {
            doCount();
        } else {
            doTile();
        }
    }

    private void outputButtonActionPerformed(ActionEvent e) {
        File chosenFile = chooseFile();
        outputField.setText(chosenFile.getAbsolutePath());
    }

    private void initUI() {
        this.igvTools = new IgvTools();
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

        zoomCombo.setSelectedIndex(IgvTools.MAX_ZOOM);
        windowSizeField.setText(String.valueOf(IgvTools.WINDOW_SIZE));
        redirectSystemStreams();

        if (mode == Mode.COVERAGE) {
            setTitle("Compute Coverage");
            label1.setText(""); //Compute coverage for an alignment or feature file");
            probeLabel.setVisible(false);
            probeField.setVisible(false);
            probeField.setVisible(false);
        } else {
            label1.setText(""); //Create a TDF file");
            setTitle("Convert to TDF");
        }
    }


    private void close() {
        System.setErr(systemErrStream);
        System.setOut(systemOutStream);
        dispose();
    }

    private void closeButtonActionPerformed(ActionEvent e) {
        System.setErr(systemErrStream);
        System.setOut(systemOutStream);
        dispose();
    }


    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        mainPanel = new JPanel();
        label1 = new JLabel();
        requiredPanel = new JPanel();
        JLabel label2 = new JLabel();
        outputLabel = new JLabel();
        outputButton = new JButton();
        inputField = new JTextField();
        inputButton = new JButton();
        outputField = new JTextField();
        genomeLabel = new JLabel();
        genomeField = new JTextField();
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
        zoomCombo = new JComboBox();
        windowSizeLabel = new JLabel();
        windowSizeField = new JTextField();
        probeLabel = new JLabel();
        probeField = new JTextField();
        probeButton = new JButton();
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

            //---- label1 ----
            label1.setText("Compute coverage for an alignment or feature file.");
            label1.setFont(new Font("Lucida Grande", Font.BOLD, 14));
            mainPanel.add(label1, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.VERTICAL,
                    new Insets(0, 0, 10, 0), 0, 0));

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
                outputButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        outputButtonActionPerformed(e);
                    }
                });
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
            }
            mainPanel.add(requiredPanel, new GridBagConstraints(1, 1, 1, 1, 1.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 10, 0), 0, 0));

            //======== tilePanel ========
            {
                tilePanel.setEnabled(true);
                tilePanel.setFont(tilePanel.getFont().deriveFont(Font.ITALIC, 10f));
                tilePanel.setBorder(new TitledBorder(null, "Advanced Options", TitledBorder.LEADING, TitledBorder.TOP));
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

                //---- zoomCombo ----
                zoomCombo.setEditable(false);
                zoomCombo.setModel(new DefaultComboBoxModel(new String[]{

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

                //---- probeLabel ----
                probeLabel.setFont(probeLabel.getFont());
                probeLabel.setToolTipText("<html>Specifies a \"bed\" file to be used to map probe identifiers to locations.  This option is useful <br>when preprocessing gct files.  The bed file should contain 4 columns: chr start end name\n<br>where name is the probe name in the gct file.");
                probeLabel.setText("Probe to Loci Mapping");
                tilePanel.add(probeLabel, new GridBagConstraints(1, 5, 1, 1, 0.0, 1.0,
                        GridBagConstraints.WEST, GridBagConstraints.NONE,
                        new Insets(0, 0, 0, 0), 0, 0));
                tilePanel.add(probeField, new GridBagConstraints(2, 5, 1, 1, 1.0, 1.0,
                        GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL,
                        new Insets(0, 0, 0, 0), 0, 0));

                //---- probeButton ----
                probeButton.setText("Browse");
                tilePanel.add(probeButton, new GridBagConstraints(3, 5, 1, 1, 0.0, 1.0,
                        GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL,
                        new Insets(0, 0, 0, 0), 0, 0));
            }
            mainPanel.add(tilePanel, new GridBagConstraints(1, 2, 1, 1, 1.0, 1.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 10, 0), 0, 0));

            //======== buttonPanel ========
            {
                buttonPanel.setLayout(new GridBagLayout());

                //---- runButton ----
                runButton.setText("Run");
                runButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        runButtonActionPerformed(e);
                    }
                });
                buttonPanel.add(runButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 1.0,
                        GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL,
                        new Insets(0, 0, 0, 0), 0, 0));

                //---- closeButton ----
                closeButton.setText("Close");
                closeButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        closeButtonActionPerformed(e);
                    }
                });
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
    private JLabel label1;
    private JPanel requiredPanel;
    private JLabel outputLabel;
    private JButton outputButton;
    private JTextField inputField;
    private JButton inputButton;
    private JTextField outputField;
    private JLabel genomeLabel;
    private JTextField genomeField;
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
    private JComboBox zoomCombo;
    private JLabel windowSizeLabel;
    private JTextField windowSizeField;
    private JLabel probeLabel;
    private JTextField probeField;
    private JButton probeButton;
    private JPanel buttonPanel;
    private JButton runButton;
    private JButton closeButton;
    private JPanel OutputPanel;
    private JScrollPane outputScroll;
    private JTextArea outputText;
    private JProgressBar progressBar;
    // JFormDesigner - End of variables declaration  //GEN-END:variables

    private void updateTextArea(final String text) {
        SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                outputText.append(text);
            }
        });
    }

    protected void redirectSystemStreams() {
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

    protected void showMessage(String tool) {
        JOptionPane.showMessageDialog(this, tool);
    }


    protected void doCount() {

        SwingWorker swingWorker = new SwingWorker() {

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

                    int extFactor = 0;
                    int strandOption = -1;

                    runButton.setEnabled(false);

                    // Check user prefs for duplicates and min mapping quality
                    int countFlags = 0;
                    boolean includeDuplicates = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_SHOW_DUPLICATES);
                    int minMappingQuality = PreferenceManager.getInstance().getAsInt(PreferenceManager.SAM_QUALITY_THRESHOLD);
                    if (includeDuplicates) {
                        countFlags += CoverageCounter.INCLUDE_DUPS;
                    }

                    int preExtFactor = 0;
                    int postExtFactor = 0;
                    igvTools.doCount(ifile, ofile, genomeId, maxZoomValue, wfs, windowSize, extFactor,
                            preExtFactor, postExtFactor, null,
                            null, minMappingQuality, countFlags);
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


    private void doTile() {

        SwingWorker swingWorker = new SwingWorker() {

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

                    String typeString = Preprocessor.getExtension("ifile");

                    runButton.setEnabled(false);
                    igvTools.toTDF(typeString, ifile, ofile, probeFile, genomeId, maxZoomValue, wfs, null, IgvTools.MAX_RECORDS_IN_RAM);
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

    protected File chooseFile() {

        //TODO Override so you can specify the file type with a string array ex: {".wig", ".tdf"}

        if (fileDialog == null) {
            fileDialog = new JFileChooser();
        }

        fileDialog.setMultiSelectionEnabled(false);
        fileDialog.setFileSelectionMode(JFileChooser.FILES_ONLY);

        int returnVal = fileDialog.showDialog(this, "Select File");
        if (returnVal == JFileChooser.CANCEL_OPTION) {
            return null;
        } else {
            File selected = fileDialog.getSelectedFile();
            return selected;
        }
    }

    public static void launch(boolean modal, String genomeId, Mode mode) {
        CoverageGui mainWindow = new CoverageGui(mode);
        mainWindow.pack();
        mainWindow.setModal(modal);
        mainWindow.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        mainWindow.setResizable(false);

        if (genomeId != null) {
            mainWindow.genomeField.setText(genomeId);
            mainWindow.genomeField.setEnabled(false);
            mainWindow.genomeField.setToolTipText("<html>To change the genome id close this window and <br>use the pulldown on the IGV batch screen.");
        }

        mainWindow.setVisible(true);
    }
}
