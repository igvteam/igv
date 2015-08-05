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
 * Created by JFormDesigner on Tue Jan 14 11:35:23 EST 2014
 */

package org.broad.igv.cursor;

import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.igv.util.encode.EncodeFileBrowser;
import org.broad.igv.util.encode.EncodeFileRecord;

import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.*;
import java.util.List;
import javax.swing.*;

/**
 * @author Stan Diamond
 */
public class CursorMainWindow extends JFrame {

    CursorModel cursorModel;  // The model

    public static void main(String[] args) {
        (new CursorMainWindow()).setVisible(true);
    }


    public CursorMainWindow() {
        initComponents();
        cursorModel = new CursorModel(this);
        cursorMainPanel1.setModel(cursorModel);
        frameWidthField.setText(String.valueOf(cursorModel.getFramePixelWidth()));
        regionSizeTextField.setText(String.valueOf(cursorModel.getFrameBPWidth()));
        ToolTipManager.sharedInstance().setDismissDelay(20000);
    }


    void updateRegionsLabel() {
        int visibleRegionCount = (int) (getWidth() / cursorModel.getFramePixelWidth()) + 1;
        final List<CursorRegion> filteredRegions = cursorModel.getFilteredRegions();
        if (filteredRegions != null) {
            regionsLabel.setText(" (" + visibleRegionCount + " / " + filteredRegions.size() + ")");
        }

    }


    private void frameWidthFieldActionPerformed(ActionEvent e) {
        try {
            double newWidth = Double.parseDouble(frameWidthField.getText().trim());
            final int regionCount = cursorModel.getFilteredRegions().size();
            if (regionCount > 0) {
                double minWidth = ((double) cursorMainPanel1.getDataPanelWidth()) / regionCount;
                newWidth = Math.max(newWidth, minWidth);
                frameWidthField.setText(String.valueOf(newWidth));
            }
            if (newWidth > 0) cursorModel.setFramePixelWidth(newWidth);
            cursorMainPanel1.repaint();
            updateRegionsLabel();
        } catch (NumberFormatException e1) {
            e1.printStackTrace();
        }
    }

    private void exitMenuItemActionPerformed(ActionEvent e) {
        setVisible(false);
        System.exit(0);
    }

    private void regionSizeTextFieldActionPerformed(ActionEvent e) {
        try {
            int newWidth = Integer.parseInt(regionSizeTextField.getText().trim());
            if (newWidth > 0) cursorModel.setFrameBPWidth(newWidth);
            cursorMainPanel1.repaint();
            updateRegionsLabel();
        } catch (NumberFormatException e1) {
            e1.printStackTrace();
        }
    }

    private void loadTracks(final Collection<EncodeFileRecord> records) {

        Runnable runnable = new Runnable() {
            public void run() {
                try {
                    startWaitCursor();
                    for (EncodeFileRecord record : records) {

                        String path = record.getPath();
                        String name = record.getTrackName();
                        Color color = null;
                        final String antibody = record.getAttributeValue("antibody");
                        if (antibody != null) {
                            color = colors.get(antibody.toUpperCase());
                        }

                        String pathLC = path.toLowerCase();
                        if (pathLC.endsWith(".gz")) pathLC = pathLC.substring(0, pathLC.length() - 3);
                        boolean loadable = pathLC.endsWith(".bed") || pathLC.endsWith(".narrowpeak") || pathLC.endsWith("broadpeak");

                        if (loadable) {
                            try {
                                CursorTrack t = CursorUtils.loadTrack(path);
                                if (t != null) {
                                    if (name != null) t.setName(name);
                                    if (color != null) t.setColor(color);
                                    cursorModel.addTrack(t);
                                    if (cursorModel.getFilteredRegions() == null || cursorModel.getFilteredRegions().isEmpty()) {
                                        cursorModel.setRegions(CursorUtils.createRegions(t));
                                    }

                                    cursorMainPanel1.addTrack(t);
                                }
                            } catch (IOException e1) {
                                e1.printStackTrace();
                            }
                        }
                    }

                    SwingUtilities.invokeLater(new Runnable() {
                        @Override
                        public void run() {
                            cursorMainPanel1.tracksAdded();
                            cursorMainPanel1.revalidate();
                            cursorMainPanel1.repaint();
                            updateRegionsLabel();
                        }
                    });

                } finally {
                    stopWaitCursor();
                }
            }
        };
        (new Thread(runnable)).start();

    }


    private void loadEncodeMenuItemActionPerformed(ActionEvent e) {

        try {
            EncodeFileBrowser browser = EncodeFileBrowser.getInstance("hg19");

            browser.setVisible(true);
            if (browser.isCanceled()) return;

            final java.util.List<EncodeFileRecord> records = browser.getSelectedRecords();
            if (records.size() > 0) {
                loadTracks(records);

            }


        } catch (IOException ex) {
            //log.error("Error opening Encode browser", e);
        }

    }


    private void startWaitCursor() {
        getGlassPane().setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
        getGlassPane().addMouseListener(nullMouseAdapter);
        getGlassPane().setVisible(true);
    }

    private void stopWaitCursor() {
        getGlassPane().setCursor(Cursor.getDefaultCursor());
        getGlassPane().removeMouseListener(nullMouseAdapter);
        getGlassPane().setVisible(false);
    }

    private void loadFIleMenuItemActionPerformed(ActionEvent e) {
        File lastDirectoryFile = PreferenceManager.getInstance().getLastTrackDirectory();


        // Tracks.  Simulates multi-file select
        File[] trackFiles = FileDialogUtils.chooseMultiple("Select Files", lastDirectoryFile, new FilenameFilter() {

            @Override
            public boolean accept(File file, String s) {
                return true;//  return file.getName().toLowerCase().endsWith("peak") || file.getName().toLowerCase().endsWith("peak.gz");
            }
        });

        if (trackFiles == null || trackFiles.length == 0) return;
        PreferenceManager.getInstance().setLastTrackDirectory(trackFiles[0]);
        List<EncodeFileRecord> records = new ArrayList<EncodeFileRecord>();
        for (File f : trackFiles) {
            records.add(new EncodeFileRecord(f.getPath(), new HashMap()));
        }
        loadTracks(records);
    }

    private void filterMenuItemActionPerformed(ActionEvent e) {

        CursorFilterDialog dlg = new CursorFilterDialog(this, cursorModel.getTracks(), cursorModel.getFilter());
        dlg.setVisible(true);
        if (!dlg.isCanceled()) {

        }
    }


    private static MouseAdapter nullMouseAdapter = new MouseAdapter() {
    };

    private static Map<String, Color> colors = new HashMap<String, Color>();

    static {
        colors = new HashMap<String, Color>();
        colors.put("H3K27AC", new Color(200, 0, 0));
        colors.put("H3K27ME3", new Color(200, 0, 0));
        colors.put("H3K36ME3", new Color(0, 0, 150));
        colors.put("H3K4ME1", new Color(0, 150, 0));
        colors.put("H3K4ME2", new Color(0, 150, 0));
        colors.put("H3K4ME3", new Color(0, 150, 0));
        colors.put("H3K9AC", new Color(100, 0, 0));
        colors.put("H3K9ME1", new Color(100, 0, 0));
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        menuBar1 = new JMenuBar();
        fileMenu = new JMenu();
        loadFileMenuItem = new JMenuItem();
        loadEncodeMenuItem = new JMenuItem();
        exitMenuItem = new JMenuItem();
        regionsMenu = new JMenu();
        filterMenuItem = new JMenuItem();
        cursorMainPanel1 = new CursorMainPanel();
        panel1 = new JPanel();
        panel2 = new JPanel();
        label2 = new JLabel();
        regionSizeTextField = new JTextField();
        panel3 = new JPanel();
        label1 = new JLabel();
        frameWidthField = new JTextField();
        regionsLabel = new JLabel();

        //======== this ========
        setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout(0, 4));

        //======== menuBar1 ========
        {

            //======== fileMenu ========
            {
                fileMenu.setText("File");

                //---- loadFileMenuItem ----
                loadFileMenuItem.setText("Load from file...");
                loadFileMenuItem.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        loadFIleMenuItemActionPerformed(e);
                    }
                });
                fileMenu.add(loadFileMenuItem);

                //---- loadEncodeMenuItem ----
                loadEncodeMenuItem.setText("Load from ENCODE...");
                loadEncodeMenuItem.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        loadEncodeMenuItemActionPerformed(e);
                    }
                });
                fileMenu.add(loadEncodeMenuItem);
                fileMenu.addSeparator();

                //---- exitMenuItem ----
                exitMenuItem.setText("Exit");
                exitMenuItem.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        exitMenuItemActionPerformed(e);
                    }
                });
                fileMenu.add(exitMenuItem);
            }
            menuBar1.add(fileMenu);

            //======== regionsMenu ========
            {
                regionsMenu.setText("Regions");

                //---- filterMenuItem ----
                filterMenuItem.setText("Filter...");
                filterMenuItem.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        filterMenuItemActionPerformed(e);
                    }
                });
                regionsMenu.add(filterMenuItem);
            }
            menuBar1.add(regionsMenu);
        }
        setJMenuBar(menuBar1);

        //---- cursorMainPanel1 ----
        cursorMainPanel1.setPreferredSize(new Dimension(135, 50));
        contentPane.add(cursorMainPanel1, BorderLayout.CENTER);

        //======== panel1 ========
        {
            panel1.setLayout(new FlowLayout(FlowLayout.RIGHT, 10, 0));

            //======== panel2 ========
            {
                panel2.setLayout(new FlowLayout(FlowLayout.CENTER, 5, 0));

                //---- label2 ----
                label2.setText("Region size (bp):");
                panel2.add(label2);

                //---- regionSizeTextField ----
                regionSizeTextField.setPreferredSize(new Dimension(60, 28));
                regionSizeTextField.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        regionSizeTextFieldActionPerformed(e);
                    }
                });
                panel2.add(regionSizeTextField);
            }
            panel1.add(panel2);

            //======== panel3 ========
            {
                panel3.setLayout(new FlowLayout(FlowLayout.CENTER, 5, 0));

                //---- label1 ----
                label1.setText("Frame width (pixels):");
                panel3.add(label1);

                //---- frameWidthField ----
                frameWidthField.setMinimumSize(new Dimension(50, 50));
                frameWidthField.setPreferredSize(new Dimension(60, 28));
                frameWidthField.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        frameWidthFieldActionPerformed(e);
                    }
                });
                panel3.add(frameWidthField);

                //---- regionsLabel ----
                regionsLabel.setHorizontalAlignment(SwingConstants.LEFT);
                regionsLabel.setMaximumSize(new Dimension(200, 0));
                regionsLabel.setPreferredSize(new Dimension(200, 28));
                regionsLabel.setHorizontalTextPosition(SwingConstants.LEFT);
                panel3.add(regionsLabel);
            }
            panel1.add(panel3);
        }
        contentPane.add(panel1, BorderLayout.SOUTH);
        setSize(1075, 770);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JMenuBar menuBar1;
    private JMenu fileMenu;
    private JMenuItem loadFileMenuItem;
    private JMenuItem loadEncodeMenuItem;
    private JMenuItem exitMenuItem;
    private JMenu regionsMenu;
    private JMenuItem filterMenuItem;
    private CursorMainPanel cursorMainPanel1;
    private JPanel panel1;
    private JPanel panel2;
    private JLabel label2;
    private JTextField regionSizeTextField;
    private JPanel panel3;
    private JLabel label1;
    private JTextField frameWidthField;
    private JLabel regionsLabel;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
