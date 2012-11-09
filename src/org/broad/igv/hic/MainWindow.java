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
 * Created by JFormDesigner on Mon Aug 02 22:04:22 EDT 2010
 */

package org.broad.igv.hic;

import com.jidesoft.swing.JideButton;
import org.apache.commons.math.linear.InvalidMatrixException;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.hic.data.*;
import org.broad.igv.hic.matrix.BasicMatrix;
import org.broad.igv.hic.matrix.DiskResidentBlockMatrix;
import org.broad.igv.hic.tools.Preprocessor;
import org.broad.igv.hic.track.EigenvectorTrack;
import org.broad.igv.hic.track.HiCTrackManager;
import org.broad.igv.hic.track.TrackPanel;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.igv.ui.util.IconFactory;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ParsingUtils;
import slider.RangeSlider;

import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.border.LineBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.dnd.*;
import java.awt.event.*;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * @author James Robinson
 */
public class MainWindow extends JFrame {

    //private static String DEFAULT_LOAD_MENU = "http://www.broadinstitute.org/igvdata/hic/hicExternalMenu.properties";
    private static String DEFAULT_LOAD_MENU = "http://iwww.broadinstitute.org/igvdata/hic/files/hicInternalMenu.properties";
    public static Color RULER_LINE_COLOR = new Color(0, 0, 0, 100);

    private ExecutorService threadExecutor = Executors.newFixedThreadPool(1);
    // The "model" object containing the state for this instance.
    private HiC hic;

    private EigenvectorTrack eigenvectorTrack;

    public static Cursor fistCursor;

    public static final int BIN_PIXEL_WIDTH = 1;

    private static MainWindow theInstance;


    public static void main(String[] args) throws IOException {

        theInstance = getInstance();
        theInstance.setVisible(true);
        //theInstance.setSize(950, 700);
        theInstance.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    }

    public static synchronized MainWindow getInstance() {
        if (theInstance == null) {
            try {
                theInstance = createMainWindow();
            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }
        return theInstance;
    }

    private MainWindow() throws IOException {

        hic = new HiC(this);

        initComponents();
        createCursors();
        pack();

        DropTarget target = new DropTarget(this, new FileDropTargetListener());
        setDropTarget(target);

        colorRangeSlider.setUpperValue(1200);
    }

    public static MainWindow createMainWindow() throws IOException {
        return new MainWindow();
    }

    public void updateColorSlider(int min, int max, int value) {
        colorRangeSlider.setMinimum(min);
        colorRangeSlider.setMaximum(max);
        colorRangeSlider.setLowerValue(0);
        colorRangeSlider.setUpperValue(value);
    }

    public void createCursors() {
        BufferedImage handImage = new BufferedImage(32, 32, BufferedImage.TYPE_INT_ARGB);

        // Make backgroun transparent
        Graphics2D g = handImage.createGraphics();
        g.setComposite(AlphaComposite.getInstance(AlphaComposite.CLEAR, 0.0f));
        Rectangle2D.Double rect = new Rectangle2D.Double(0, 0, 32, 32);
        g.fill(rect);

        // Draw hand image in middle
        g = handImage.createGraphics();
        g.drawImage(IconFactory.getInstance().getIcon(IconFactory.IconID.FIST).getImage(), 0, 0, null);
        MainWindow.fistCursor = getToolkit().createCustomCursor(handImage, new Point(8, 6), "Move");
    }

    public HeatmapPanel getHeatmapPanel() {
        return heatmapPanel;
    }


    public void updateZoom(HiC.Unit unit, int newZoom) {

        int idx = unit == HiC.Unit.BP ? newZoom : newZoom + 9;
        resolutionSlider.setValue(idx);
    }

    /**
     * Chromosome "0" is whole genome
     *
     * @param chromosomes
     */
    public void setChromosomes(Chromosome[] chromosomes) {
        hic.setChromosomes(chromosomes);
        int[] chromosomeBoundaries = new int[chromosomes.length - 1];
        long bound = 0;
        for (int i = 1; i < chromosomes.length; i++) {
            Chromosome c = chromosomes[i];
            bound += (c.getLength() / 1000);
            chromosomeBoundaries[i - 1] = (int) bound;
        }
        heatmapPanel.setChromosomeBoundaries(chromosomeBoundaries);
    }


    public void setSelectedChromosomes(Chromosome xChrom, Chromosome yChrom) {
        chrBox1.setSelectedIndex(yChrom.getIndex());
        chrBox2.setSelectedIndex(xChrom.getIndex());
        refreshChromosomes();
    }


    private void load(String file) throws IOException {
        if (file.endsWith("hic")) {

            String magicString = DatasetReaderV2.getMagicString(file);

            DatasetReader reader;
            if (magicString.equals("HIC")) {
                reader = new DatasetReaderV2(file);
            } else {
                reader = new DatasetReaderV1(file);
                // file not actually read, usually canceled the read of password-protected file
                if (reader.getVersion() == -1)
                    return;
            }

            hic.dataset = reader.read();
            if (hic.dataset.getVersion() <= 1) {
                JOptionPane.showMessageDialog(this, "This version of \"hic\" format is no longer supported");
                return;
            }

            setChromosomes(hic.dataset.getChromosomes());
            chrBox1.setModel(new DefaultComboBoxModel(hic.getChromosomes()));
            chrBox2.setModel(new DefaultComboBoxModel(hic.getChromosomes()));

            displayOptionComboBox.setModel(new DefaultComboBoxModel(new DisplayOption[]{
                    DisplayOption.OBSERVED,
                    DisplayOption.OE,
                    DisplayOption.PEARSON}));

            displayOptionComboBox.setSelectedIndex(0);
            //setTitle(file);
            hic.xContext = null;
            hic.yContext = null;
            refreshChromosomes();
        } else {
            // error -- unknown file type
        }

    }

    private void refreshChromosomes() {


        Runnable runnable = new Runnable() {
            public void run() {
                Chromosome chr1 = (Chromosome) chrBox1.getSelectedItem();
                Chromosome chr2 = (Chromosome) chrBox2.getSelectedItem();

                int t1 = chr1.getIndex();
                int t2 = chr2.getIndex();

                if (t1 > t2) {
                    Chromosome tmp = chr2;
                    chr2 = chr1;
                    chr1 = tmp;
                }

                if (chr2.getName().equals("All") ||
                        chr2 != hic.xContext.getChromosome() || chr1 != hic.yContext.getChromosome()) {

                    hic.xContext = new Context(chr2);
                    hic.yContext = new Context(chr1);
                    rulerPanel2.setFrame(hic.xContext, HiCRulerPanel.Orientation.HORIZONTAL);
                    rulerPanel1.setFrame(hic.yContext, HiCRulerPanel.Orientation.VERTICAL);

                    hic.matrix = hic.dataset.getMatrix(chr1, chr2);
                    if (hic.matrix == null) {
                        //?? TODO -- this is probably an error
                    } else {
                        setInitialZoom();  // TODO -- Will probably trigger a repaint
                    }
                }
                if (t1 != t2 || chr2.getName().equals("All")) { // || hic.getZoomToDensityMap() == null) {
                    viewEigenvector.setEnabled(false);
                    displayOptionComboBox.setSelectedIndex(0);
                    displayOptionComboBox.setEnabled(false);
                } else {
                    viewEigenvector.setEnabled(true);
                    displayOptionComboBox.setEnabled(true);
                }

                refresh();
            }
        };

        executeLongRunningTask(runnable);

    }

    void refresh() {
        getHeatmapPanel().clearTileCache();
        repaint();
        if (hic.matrix != null) {
            //MatrixZoomData.ScaleParameters scaleParameters = zd.computeScaleParameters();
            // Refresh the thumbnail, using lowest resolution data
            // TODO -- this only needs to be done in some circumstances

            MatrixZoomData zd0 = hic.matrix.getObservedMatrix(0);
            Image thumbnail = heatmapPanel.getThumbnailImage(
                    zd0,
                    thumbnailPanel.getWidth(),
                    thumbnailPanel.getHeight(),
                    hic.getDisplayOption());
            thumbnailPanel.setImage(thumbnail);
        }
    }

    private void setInitialZoom() {

        if (hic.xContext.getChromosome().getName().equals("All")) {
            resolutionSlider.setValue(0);
            resolutionSlider.setEnabled(false); //
            hic.setZoom(0, -1, -1);
        } else {
            resolutionSlider.setEnabled(true);
            int initialZoom = 1;

//            Find right zoom level
//            int pixels = getHeatmapPanel().getWidth();
//            int len = (Math.max(xContext.getChrLength(), yContext.getChrLength()));
//            int maxNBins = pixels / BIN_PIXEL_WIDTH;
//            int bp_bin = len / maxNBins;
//            int initalZoom = HiCGlobals.zoomBinSizes.length - 1;
//            for (int z = 1; z < HiCGlobals.zoomBinSizes.length; z++) {
//                if (HiCGlobals.zoomBinSizes[z] < bp_bin) {
//                    initialZoom = z - 1;
//                    break;
//                }
//            }
            resolutionSlider.setValue(initialZoom);
            hic.setZoom(initialZoom, -1, -1);
        }
    }

    private void refreshButtonActionPerformed(ActionEvent e) {
        refreshChromosomes();
    }

    private void loadMenuItemActionPerformed(ActionEvent e) {
        FileDialog dlg = new FileDialog(this);
        dlg.setMode(FileDialog.LOAD);
        dlg.setVisible(true);
        String file = dlg.getFile();
        if (file != null) {
            try {
                File f = new File(dlg.getDirectory(), dlg.getFile());
                load(f.getAbsolutePath());
            } catch (IOException e1) {
                e1.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }
    }


    private void loadFromURLActionPerformed(ActionEvent e) {
        String url = JOptionPane.showInputDialog("Enter URL: ");
        if (url != null) {
            try {
                load(url);
            } catch (IOException e1) {
                e1.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }
    }

    private void exitActionPerformed(ActionEvent e) {
        setVisible(false);
        dispose();
        System.exit(0);
    }

    private void colorRangeSliderStateChanged(ChangeEvent e) {
        int min = colorRangeSlider.getLowerValue();
        int max = colorRangeSlider.getUpperValue();
        heatmapPanel.setObservedRange(min, max);
    }

    private void chrBox1ActionPerformed(ActionEvent e) {
        if (chrBox1.getSelectedIndex() == 0) {
            chrBox2.setSelectedIndex(0);
        }
    }

    private void chrBox2ActionPerformed(ActionEvent e) {
        if (chrBox2.getSelectedIndex() == 0) {
            chrBox1.setSelectedIndex(0);
        }
    }

    private void displayOptionComboBoxActionPerformed(ActionEvent e) {

        DisplayOption option = (DisplayOption) (displayOptionComboBox.getSelectedItem());
        hic.setDisplayOption(option);
        switch (option) {
            case OBSERVED:
                break;
            case OE:
                break;
            case PEARSON:
                BasicMatrix bm = hic.zd.getPearsons();
                if (bm != null) {
                    float lv = bm.getLowerValue();
                    float uv = bm.getUpperValue();

                    // colorRangeSlider.setLowerValue(lv);
                    // colorRangeSlider.setUpperValue(uv);
                }

        }


    }


    void updateEigenvectorTrack() {
        boolean show = viewEigenvector.isSelected();
        if (show) {
            trackPanel.setEigenvectorTrack(eigenvectorTrack);
            eigenvectorTrack.forceRefresh();
        }
        updateTrackPanel();
    }

    public void setViewEigenvector(boolean flag) {
        viewEigenvector.setSelected(flag);
    }


    private void getEigenvectorActionPerformed(ActionEvent e) {
        double[] rv;
        try {
            String number = JOptionPane.showInputDialog("Which eigenvector do you want to see?");
            int num = Integer.parseInt(number) - 1;

            rv = hic.getEigenvector(num);

            if (rv != null) {
                String str = "";
                for (int i = 0; i < rv.length; i++) {
                    str += rv[i] + "\n";

                }
                JTextArea textArea = new JTextArea(str, 20, 20);
                textArea.setEditable(false);
                textArea.selectAll();
                JScrollPane pane = new JScrollPane(textArea);
                JFrame frame = new JFrame("Principal Eigenvector");
                frame.getContentPane().add(pane);
                frame.pack();
                frame.setVisible(true);

            } else {
                System.err.println("No densities available for this file.");
            }

        } catch (InvalidMatrixException error) {
            JOptionPane.showMessageDialog(this, "Unable to calculate eigenvectors after 30 iterations",
                    "Eigenvector error", JOptionPane.ERROR_MESSAGE);
        } catch (NumberFormatException error) {
            JOptionPane.showMessageDialog(this, "You must enter a valid number.\n" + error.getMessage(),
                    "Eigenvector error", JOptionPane.ERROR_MESSAGE);
        }


    }


    /**
     * Untility function to execute a task in a worker thread.  The method is on MainWindow because the glassPane
     * is used to display a wait cursor and block events.
     *
     * @param runnable
     * @return
     */
    public Future executeLongRunningTask(final Runnable runnable) {

        Callable<Object> wrapper = new Callable<Object>() {
            public Object call() throws Exception {
                final Component glassPane = showGlassPane();
                try {
                    runnable.run();
                    return "done";
                } finally {
                    glassPane.setVisible(false);
                }

            }
        };

        return threadExecutor.submit(wrapper);
    }

    public Component showGlassPane() {
        final Component glassPane = getGlassPane();
        glassPane.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
        glassPane.setVisible(true);
        return glassPane;
    }

    public void hideGlassPane() {
        getGlassPane().setVisible(false);
    }

    public void updateTrackPanel() {
        if (HiCTrackManager.getLoadedTracks().size() > 0 || viewEigenvector.isSelected()) {
            trackPanel.setVisible(true);
        } else {
            trackPanel.setVisible(false);
        }
        invalidate();
        pack();
        repaint();

    }

    /**
     * Listener for drag&drop actions
     */
    class FileDropTargetListener implements DropTargetListener {


        public FileDropTargetListener() {
        }

        public void dragEnter(DropTargetDragEvent event) {

            if (!isDragAcceptable(event)) {
                event.rejectDrag();
                return;
            }
        }

        public void dragExit(DropTargetEvent event) {
        }

        public void dragOver(DropTargetDragEvent event) {
            // you can provide visual feedback here
        }

        public void dropActionChanged(DropTargetDragEvent event) {
            if (!isDragAcceptable(event)) {
                event.rejectDrag();
                return;
            }
        }

        public void drop(DropTargetDropEvent event) {
            if (!isDropAcceptable(event)) {
                event.rejectDrop();
                return;
            }

            event.acceptDrop(DnDConstants.ACTION_COPY);

            Transferable transferable = event.getTransferable();

            try {
                java.util.List<File> files = (java.util.List<File>) transferable.getTransferData(DataFlavor.javaFileListFlavor);
                load(files.get(0).getAbsolutePath());

            } catch (Exception e) {
                String obj = null;
                try {
                    obj = transferable.getTransferData(DataFlavor.stringFlavor).toString();
                    if (HttpUtils.isRemoteURL(obj)) {
                        load(obj);
                    }
                } catch (Exception e1) {
                }

            }
            repaint();
            event.dropComplete(true);
        }


        public boolean isDragAcceptable(DropTargetDragEvent event) {
            //  Check the  available data flavors here
            //  Currently accepting all flavors
            return (event.getDropAction() & DnDConstants.ACTION_COPY_OR_MOVE) != 0;
        }

        public boolean isDropAcceptable(DropTargetDropEvent event) {
            //  Check the  available data flavors here
            //  Currently accepting all flavors
            return (event.getDropAction() & DnDConstants.ACTION_COPY_OR_MOVE) != 0;
        }
    }

    private boolean addPredefinedLoadItems(JMenu fileMenu) {
        InputStream is = null;
        Properties properties = null;
        boolean isInternal = false;

        try {
            String url = System.getProperty("loadMenu");
            if (url == null) url = DEFAULT_LOAD_MENU;
            is = ParsingUtils.openInputStream(url);
            properties = new Properties();
            properties.load(is);
            isInternal = (url == DEFAULT_LOAD_MENU);
        } catch (Exception error) {
            System.err.println("Can't find mainwindow.properties.");
            return false;
        }
        // TreeSet is sorted, so properties file is implemented in order
        TreeSet<String> keys = new TreeSet<String>(properties.stringPropertyNames());
        for (String key : keys) {
            String value = properties.getProperty(key);
            final String[] values = value.split(",");
            if (values.length != 3 && values.length != 1) {
                System.err.println("Improperly formatted mainwindow.properties file");
                return false;
            }
            if (values.length == 1) {
                fileMenu.addSeparator();
            } else {
                final int maxValue = Integer.parseInt(values[2]);
                final String name = values[0];
                JMenuItem item = new JMenuItem(values[0]);
                item.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        try {
                            heatmapPanel.setObservedRange(0, maxValue);
                            colorRangeSlider.setMaximum(maxValue);
                            colorRangeSlider.setMinimum(0);
                            //colorRangeSlider.setMajorTickSpacing((int) (maxValue / 10));
                            colorRangeSlider.setUpperValue((int) (maxValue * 3 / 4));
                            hic.reset();
                            load(values[1]);
                            setTitle(name);
                        } catch (IOException e1) {
                            JOptionPane.showMessageDialog(MainWindow.this, "Error loading data: " + e1.getMessage());
                        }

                    }
                });
                fileMenu.add(item);
            }
        }
        fileMenu.addSeparator();
        return isInternal;
    }


    private void initComponents() {
        JPanel mainPanel = new JPanel();
        JPanel toolbarPanel = new JPanel();


        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());
        mainPanel.setLayout(new BorderLayout());

        toolbarPanel.setBorder(null);
        toolbarPanel.setLayout(new GridLayout());


        // --- Chromosome panel ---
        JPanel chrSelectionPanel = new JPanel();
        chrSelectionPanel.setBorder(LineBorder.createGrayLineBorder());
        chrSelectionPanel.setMinimumSize(new Dimension(130, 57));
        chrSelectionPanel.setPreferredSize(new Dimension(130, 57));
        chrSelectionPanel.setLayout(new BorderLayout());

        JPanel chrLabelPanel = new JPanel();
        JLabel chrLabel = new JLabel("Chromosomes");
        chrLabel.setHorizontalAlignment(SwingConstants.CENTER);
        chrLabelPanel.setBackground(new Color(204, 204, 204));
        chrLabelPanel.setLayout(new BorderLayout());
        chrLabelPanel.add(chrLabel, BorderLayout.CENTER);
        chrSelectionPanel.add(chrLabelPanel, BorderLayout.PAGE_START);

        JPanel chrButtonPanel = new JPanel();
        chrButtonPanel.setBackground(new Color(238, 238, 238));
        chrButtonPanel.setLayout(new BoxLayout(chrButtonPanel, BoxLayout.X_AXIS));

        //---- chrBox1 ----
        chrBox1 = new JComboBox();
        chrBox1.setModel(new DefaultComboBoxModel(new String[]{"All"}));
        chrBox1.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                chrBox1ActionPerformed(e);
            }
        });
        chrButtonPanel.add(chrBox1);

        //---- chrBox2 ----
        chrBox2 = new JComboBox();
        chrBox2.setModel(new DefaultComboBoxModel(new String[]{"All"}));
        chrBox2.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                chrBox2ActionPerformed(e);
            }
        });
        chrButtonPanel.add(chrBox2);

        //---- refreshButton ----
        JideButton refreshButton = new JideButton();
        refreshButton.setIcon(new ImageIcon(getClass().getResource("/toolbarButtonGraphics/general/Refresh24.gif")));
        refreshButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                refreshButtonActionPerformed(e);
            }
        });
        chrButtonPanel.add(refreshButton);

        chrSelectionPanel.add(chrButtonPanel, BorderLayout.CENTER);

        toolbarPanel.add(chrSelectionPanel);

        //======== displayOptionPanel ========
        JPanel displayOptionPanel = new JPanel();

        displayOptionPanel.setBackground(new Color(238, 238, 238));
        displayOptionPanel.setBorder(LineBorder.createGrayLineBorder());
        displayOptionPanel.setLayout(new BorderLayout());

        JPanel displayOptionLabelPanel = new JPanel();
        displayOptionLabelPanel.setBackground(new Color(204, 204, 204));
        displayOptionLabelPanel.setLayout(new BorderLayout());

        JLabel displayOptionLabel = new JLabel("Show");
        displayOptionLabel.setHorizontalAlignment(SwingConstants.CENTER);
        displayOptionLabelPanel.add(displayOptionLabel, BorderLayout.CENTER);

        displayOptionPanel.add(displayOptionLabelPanel, BorderLayout.PAGE_START);

        JPanel displayOptionButtonPanel = new JPanel();
        displayOptionButtonPanel.setBorder(new EmptyBorder(0, 10, 0, 10));
        displayOptionButtonPanel.setLayout(new GridLayout(1, 0, 20, 0));
        displayOptionComboBox = new JComboBox();
        displayOptionComboBox.setModel(new DefaultComboBoxModel(new String[]{DisplayOption.OBSERVED.toString()}));
        displayOptionComboBox.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                displayOptionComboBoxActionPerformed(e);
            }
        });
        displayOptionButtonPanel.add(displayOptionComboBox);

        displayOptionPanel.add(displayOptionButtonPanel, BorderLayout.CENTER);

        toolbarPanel.add(displayOptionPanel);

        //======== colorRangePanel ========

        JPanel colorRangePanel = new JPanel();
        colorRangeSlider = new RangeSlider();
        colorRangePanel.setBorder(LineBorder.createGrayLineBorder());
        colorRangePanel.setMinimumSize(new Dimension(96, 70));
        colorRangePanel.setPreferredSize(new Dimension(202, 70));
        colorRangePanel.setMaximumSize(new Dimension(32769, 70));
        colorRangePanel.setLayout(new BorderLayout());


        //---- colorRangeLabel ----
        JLabel colorRangeLabel = new JLabel("Color Range");
        colorRangeLabel.setHorizontalAlignment(SwingConstants.CENTER);
        colorRangeLabel.setToolTipText("Range of color scale in counts per mega-base squared.");
        colorRangeLabel.setHorizontalTextPosition(SwingConstants.CENTER);
        colorRangeLabel.addMouseListener(new MouseAdapter() {
            @Override
            public void mousePressed(MouseEvent e) {
                if (e.isPopupTrigger()) {
                    ColorRangeDialog rangeDialog = new ColorRangeDialog(MainWindow.this, colorRangeSlider);
                    rangeDialog.setVisible(true);
                }
            }

            @Override
            public void mouseClicked(MouseEvent e) {
                ColorRangeDialog rangeDialog = new ColorRangeDialog(MainWindow.this, colorRangeSlider);
                rangeDialog.setVisible(true);
            }
        });
        JPanel colorLabelPanel = new JPanel();
        colorLabelPanel.setBackground(new Color(204, 204, 204));
        colorLabelPanel.setLayout(new BorderLayout());
        colorLabelPanel.add(colorRangeLabel, BorderLayout.CENTER);

        colorRangePanel.add(colorLabelPanel, BorderLayout.PAGE_START);

        //---- colorRangeSlider ----
        colorRangeSlider.setPaintTicks(false);
        colorRangeSlider.setPaintLabels(false);
        //colorRangeSlider.setMinorTickSpacing(100);
        colorRangeSlider.setMaximumSize(new Dimension(32767, 52));
        colorRangeSlider.setPreferredSize(new Dimension(200, 52));
        colorRangeSlider.setMinimumSize(new Dimension(36, 52));

        colorRangeSlider.setMaximum(2000);
        colorRangeSlider.setLowerValue(0);
        //colorRangeSlider.setPaintTicks(true);
        //colorRangeSlider.setMajorTickSpacing(2000);
        colorRangeSlider.setUpperValue(500);
        colorRangeSlider.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                colorRangeSliderStateChanged(e);
            }
        });
        colorRangePanel.add(colorRangeSlider, BorderLayout.PAGE_END);

//        JPanel colorRangeTextPanel = new JPanel();
//        colorRangeTextPanel.setLayout(new FlowLayout());
//        JTextField minField = new JTextField();
//        minField.setPreferredSize(new Dimension(50, 15));
//        colorRangeTextPanel.add(minField);
//        colorRangeTextPanel.add(new JLabel(" - "));
//        JTextField maxField = new JTextField();
//        maxField.setPreferredSize(new Dimension(50, 15));
//        colorRangeTextPanel.add(maxField);
//        colorRangeTextPanel.setPreferredSize(new Dimension(200, 52));
//        colorRangePanel.add(colorRangeTextPanel, BorderLayout.PAGE_END);


        toolbarPanel.add(colorRangePanel);

        //======== resolutionPanel ========

        JPanel resolutionPanel = new JPanel();
        resolutionPanel.setBorder(LineBorder.createGrayLineBorder());
        resolutionPanel.setLayout(new BorderLayout());

        JLabel resolutionLabel = new JLabel("Resolution");
        resolutionLabel.setHorizontalAlignment(SwingConstants.CENTER);
        resolutionLabel.setBackground(new Color(204, 204, 204));

        JPanel resolutionLabelPanel = new JPanel();
        resolutionLabelPanel.setBackground(new Color(204, 204, 204));
        resolutionLabelPanel.setLayout(new BorderLayout());
        resolutionLabelPanel.add(resolutionLabel, BorderLayout.CENTER);

        resolutionPanel.add(resolutionLabelPanel, BorderLayout.PAGE_START);


        JPanel resolutionButtonPanel = new JPanel();
        resolutionButtonPanel.setLayout(new BoxLayout(resolutionButtonPanel, BoxLayout.X_AXIS));

        //---- resolutionSlider ----
        resolutionSlider = new JSlider();
        resolutionSlider.setMaximum(9);
        resolutionSlider.setMajorTickSpacing(1);
        resolutionSlider.setPaintTicks(true);
        resolutionSlider.setSnapToTicks(true);
        resolutionSlider.setPaintLabels(true);
        resolutionSlider.setMinorTickSpacing(1);


        // TODO -- the available resolutions should be read from the dataset (hic) file
        Dictionary<Integer, JLabel> resolutionLabels = new Hashtable<Integer, JLabel>();
        Font f = FontManager.getFont(8);
        for (int i = 0; i < Preprocessor.bpResLabels.length; i++) {
            if ((i + 1) % 2 == 0) {
                final JLabel tickLabel = new JLabel(Preprocessor.bpResLabels[i]);
                tickLabel.setFont(f);
                resolutionLabels.put(i, tickLabel);
            }
        }

        JLabel fragLabel = new JLabel("1f");
        fragLabel.setFont(f);
        resolutionLabels.put(9, fragLabel);
        resolutionSlider.setLabelTable(resolutionLabels);


        // Setting the zoom should always be done by calling resolutionSlider.setValue() so work isn't done twice.
        resolutionSlider.addChangeListener(new ChangeListener() {
            // Change zoom level while staying centered on current location.
            // Centering is relative to the bounds of the data, which might not be the bounds of the window

            public void stateChanged(ChangeEvent e) {
                if (!resolutionSlider.getValueIsAdjusting()) {


                    int idx = resolutionSlider.getValue();

                    // Temporary hack -- "9" means fragment resolution  (1f).
                    HiC.Unit unit = idx == 9 ? HiC.Unit.FRAG : HiC.Unit.BP;

                    int tmp = idx >= 9 ? idx - 9 : idx;
                    // Temporary hacks for non-internal hic viewer; called when set maximum called, also when drawing initially

                    if (hic.dataset == null) {
                        Dictionary<Integer, JLabel> resolutionLabels = new Hashtable<Integer, JLabel>();
                        Font f = FontManager.getFont(8);
                        for (int i = 0; i < 7; i++) {
                            final JLabel tickLabel = new JLabel(Preprocessor.bpResLabels[i]);
                            tickLabel.setFont(f);
                            resolutionLabels.put(i, tickLabel);
                        }
                        resolutionSlider.setLabelTable(resolutionLabels);
                        return;
                    }

                    int zoom = Math.max(0, Math.min(tmp, hic.dataset.getNumberZooms(unit)));

                    if (hic.zd != null && zoom == hic.zd.getZoom() && unit == hic.getUnit()) {
                        // Nothing to do
                        return;
                    }

                    if (hic.xContext != null) {

                        hic.setUnit(unit);
                        hic.xContext.setScaleFactor(1.0);
                        hic.yContext.setScaleFactor(1.0);

                        //int centerBinX = hic.xContext.getBinOrigin() + heatmapPanel.getWidth() / 2;
                        //int centerBinY = hic.yContext.getBinOrigin() + heatmapPanel.getHeight() / 2;
                        int centerBinX = hic.xContext.getBinOrigin() + (int) (heatmapPanel.getWidth() / (2 * hic.xContext.getScaleFactor()));
                        int centerBinY = hic.yContext.getBinOrigin() + (int) (heatmapPanel.getHeight() / (2 * hic.yContext.getScaleFactor()));

                        if (hic.zd == null) {
                            hic.setZoom(zoom, 0, 0);
                        } else {


                            int xGenome = hic.zd.getxGridAxis().getGenomicMid(centerBinX);
                            int yGenome = hic.zd.getyGridAxis().getGenomicMid(centerBinY);
                            hic.setZoom(zoom, xGenome, yGenome);
                        }
                    }
                    //zoomInButton.setEnabled(newZoom < MAX_ZOOM);
                    //zoomOutButton.setEnabled(newZoom > 0);
                }
            }
        });
        resolutionButtonPanel.add(resolutionSlider);

        resolutionPanel.add(resolutionButtonPanel, BorderLayout.CENTER);

        toolbarPanel.add(resolutionPanel);

        mainPanel.add(toolbarPanel, BorderLayout.NORTH);


        //======== hiCPanel ========


        final JPanel hiCPanel = new JPanel();
        hiCPanel.setLayout(new HiCLayout());

        //---- rulerPanel2 ----
        rulerPanel2 = new HiCRulerPanel(hic);
        rulerPanel2.setMaximumSize(new Dimension(4000, 50));
        rulerPanel2.setMinimumSize(new Dimension(1, 50));
        rulerPanel2.setPreferredSize(new Dimension(1, 50));
        rulerPanel2.setBorder(null);

        JPanel panel2_5 = new JPanel();
        panel2_5.setLayout(new BorderLayout());
        panel2_5.add(rulerPanel2, BorderLayout.SOUTH);

        trackPanel = new TrackPanel(hic);
        trackPanel.setMaximumSize(new Dimension(4000, 50));
        trackPanel.setPreferredSize(new Dimension(1, 50));
        trackPanel.setMinimumSize(new Dimension(1, 50));
        trackPanel.setBorder(LineBorder.createGrayLineBorder());
        trackPanel.setVisible(false);
        panel2_5.add(trackPanel, BorderLayout.NORTH);

        hiCPanel.add(panel2_5, BorderLayout.NORTH);


        //---- rulerPanel1 ----
        rulerPanel1 = new HiCRulerPanel(hic);
        rulerPanel1.setMaximumSize(new Dimension(50, 4000));
        rulerPanel1.setPreferredSize(new Dimension(50, 700));
        rulerPanel1.setBorder(null);
        rulerPanel1.setMinimumSize(new Dimension(50, 1));
        hiCPanel.add(rulerPanel1, BorderLayout.WEST);

        //---- heatmapPanel ----
        heatmapPanel = new HeatmapPanel(this, hic);
        heatmapPanel.setBorder(LineBorder.createBlackLineBorder());
        heatmapPanel.setMaximumSize(new Dimension(700, 700));
        heatmapPanel.setMinimumSize(new Dimension(700, 700));
        heatmapPanel.setPreferredSize(new Dimension(700, 700));
        heatmapPanel.setBackground(new Color(238, 238, 238));
        hiCPanel.add(heatmapPanel, BorderLayout.CENTER);


        //======== panel8 ========

        JPanel rightSidePanel = new JPanel();
        rightSidePanel.setMaximumSize(new Dimension(120, 100));
        rightSidePanel.setBorder(new EmptyBorder(0, 10, 0, 0));
        rightSidePanel.setLayout(null);

        //---- thumbnailPanel ----
        thumbnailPanel = new ThumbnailPanel(this, hic);
        thumbnailPanel.setMaximumSize(new Dimension(100, 100));
        thumbnailPanel.setMinimumSize(new Dimension(100, 100));
        thumbnailPanel.setPreferredSize(new Dimension(100, 100));
        thumbnailPanel.setBorder(LineBorder.createBlackLineBorder());
        thumbnailPanel.setPreferredSize(new Dimension(100, 100));
        thumbnailPanel.setBounds(new Rectangle(new Point(20, 0), thumbnailPanel.getPreferredSize()));
        rightSidePanel.add(thumbnailPanel);

        //======== xPlotPanel ========

        xPlotPanel = new JPanel();
        xPlotPanel.setPreferredSize(new Dimension(250, 100));
        xPlotPanel.setLayout(null);

        rightSidePanel.add(xPlotPanel);
        xPlotPanel.setBounds(10, 100, xPlotPanel.getPreferredSize().width, 228);

        //======== yPlotPanel ========

        yPlotPanel = new JPanel();
        yPlotPanel.setPreferredSize(new Dimension(250, 100));
        yPlotPanel.setLayout(null);

        rightSidePanel.add(yPlotPanel);
        yPlotPanel.setBounds(10, 328, yPlotPanel.getPreferredSize().width, 228);

        // compute preferred size
        Dimension preferredSize = new Dimension();
        for (int i = 0; i < rightSidePanel.getComponentCount(); i++) {
            Rectangle bounds = rightSidePanel.getComponent(i).getBounds();
            preferredSize.width = Math.max(bounds.x + bounds.width, preferredSize.width);
            preferredSize.height = Math.max(bounds.y + bounds.height, preferredSize.height);
        }
        Insets insets = rightSidePanel.getInsets();
        preferredSize.width += insets.right;
        preferredSize.height += insets.bottom;
        rightSidePanel.setMinimumSize(preferredSize);
        rightSidePanel.setPreferredSize(preferredSize);


        hiCPanel.add(rightSidePanel, BorderLayout.EAST);

        mainPanel.add(hiCPanel, BorderLayout.CENTER);

        contentPane.add(mainPanel, BorderLayout.CENTER);

        JMenuBar menuBar = createMenuBar(hiCPanel);
        contentPane.add(menuBar, BorderLayout.NORTH);

        // setup the glass pane to display a wait cursor when visible, and to grab all mouse events
        rootPane.getGlassPane().setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
        rootPane.getGlassPane().addMouseListener(new MouseAdapter() {
        });

    }


    private JMenuBar createMenuBar(final JPanel hiCPanel) {


        JMenuBar menuBar = new JMenuBar();

        //======== fileMenu ========
        JMenu fileMenu = new JMenu("File");


        //---- loadMenuItem ----
        JMenuItem loadMenuItem = new JMenuItem();
        loadMenuItem.setText("Load...");
        loadMenuItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                loadMenuItemActionPerformed(e);
            }
        });
        fileMenu.add(loadMenuItem);

        //---- loadFromURL ----
        JMenuItem loadFromURL = new JMenuItem();
        //JMenuItem getEigenvector = new JMenuItem();
        //final JCheckBoxMenuItem viewDNAseI;

        loadFromURL.setText("Load from URL ...");
        loadFromURL.setName("loadFromURL");
        loadFromURL.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                loadFromURLActionPerformed(e);
            }
        });
        fileMenu.add(loadFromURL);
        fileMenu.addSeparator();

        // Pre-defined datasets.  TODO -- generate from a file
        boolean isInternal = addPredefinedLoadItems(fileMenu);

        JMenuItem saveToImage = new JMenuItem();
        saveToImage.setText("Save to image");
        saveToImage.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                BufferedImage image = (BufferedImage) createImage(1000, 1000);
                Graphics g = image.createGraphics();
                hiCPanel.paint(g);

                JFileChooser fc = new JFileChooser();
                fc.setSelectedFile(new File("image.png"));
                int actionDialog = fc.showSaveDialog(null);
                if (actionDialog == JFileChooser.APPROVE_OPTION) {
                    File file = fc.getSelectedFile();
                    if (file.exists()) {
                        actionDialog = JOptionPane.showConfirmDialog(null, "Replace existing file?");
                        if (actionDialog == JOptionPane.NO_OPTION || actionDialog == JOptionPane.CANCEL_OPTION)
                            return;
                    }
                    try {
                        // default if they give no format or invalid format
                        String fmt = "jpg";
                        int ind = file.getName().indexOf(".");
                        if (ind != -1) {
                            String ext = file.getName().substring(ind + 1);
                            String[] strs = ImageIO.getWriterFormatNames();
                            for (String aStr : strs)
                                if (ext.equals(aStr))
                                    fmt = ext;
                        }
                        ImageIO.write(image.getSubimage(0, 0, 600, 600), fmt, file);
                    } catch (IOException ie) {
                        System.err.println("Unable to write " + file + ": " + ie);
                    }
                }
            }
        });
        fileMenu.add(saveToImage);
        //getEigenvector = new JMenuItem("Get principal eigenvector");
        //getEigenvector.addActionListener(new ActionListener() {
        //    public void actionPerformed(ActionEvent e) {
        //        getEigenvectorActionPerformed(e);
        //    }
        //});
        //fileMenu.add(getEigenvector);
        //---- exit ----
        JMenuItem exit = new JMenuItem();
        exit.setText("Exit");
        exit.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                exitActionPerformed(e);
            }
        });
        fileMenu.add(exit);


        menuBar.add(fileMenu);

        //======== Tracks menu ========

        JMenu tracksMenu = new JMenu("Tracks");

        viewEigenvector = new JCheckBoxMenuItem("View Eigenvector...");
        viewEigenvector.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                if (viewEigenvector.isSelected()) {
                    if (eigenvectorTrack == null) {
                        eigenvectorTrack = new EigenvectorTrack("eigen", "Eigenvectors", hic);
                        trackPanel.setEigenvectorTrack(eigenvectorTrack);
                    }
                } else {
                    trackPanel.setEigenvectorTrack(null);
                    if (HiCTrackManager.getLoadedTracks().isEmpty()) {
                        trackPanel.setVisible(false);
                    }
                }
                updateTrackPanel();
            }
        });
        viewEigenvector.setEnabled(false);
        tracksMenu.add(viewEigenvector);

        JMenuItem loadItem = new JMenuItem("Load...");
        loadItem.addActionListener(new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                HiCTrackManager.loadHostedTrack(MainWindow.this, hic);
            }

        });
        tracksMenu.add(loadItem);

        JMenuItem loadFromFileItem = new JMenuItem("Load from file...");
        loadFromFileItem.addActionListener(new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                HiCTrackManager.loadTrackFromFile(MainWindow.this, hic);
            }

        });
        tracksMenu.add(loadFromFileItem);

        menuBar.add(tracksMenu);

        //======== Extras menu ========
        JMenu extrasMenu = new JMenu("Extras");

        JMenuItem dumpPearsons = new JMenuItem("Dump pearsons matrix ...");
        dumpPearsons.addActionListener(new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                BasicMatrix pearsons = hic.zd.getPearsons();
                try {
                    String chr1 = hic.getChromosomes()[hic.zd.getChr1()].getName();
                    String chr2 = hic.getChromosomes()[hic.zd.getChr2()].getName();
                    int binSize = hic.zd.getBinSize();
                    File initFile = new File("pearsons_" + chr1 + "_" + "_" + chr2 + "_" + binSize + ".bin");
                    File f = FileDialogUtils.chooseFile("Save pearsons", null, initFile, FileDialogUtils.SAVE);
                    if (f != null) {
                        ScratchPad.dumpPearsonsBinary(pearsons, chr1, chr2, hic.zd.getBinSize(), f);
                    }
                } catch (IOException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }

        });


        JMenuItem dumpEigenvector = new JMenuItem("Dump eigenvector ...");
        dumpEigenvector.addActionListener(new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                try {
                    ScratchPad.dumpEigenvector(hic);
                } catch (IOException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }

        });
        extrasMenu.add(dumpEigenvector);


        JMenuItem readPearsons = new JMenuItem("Read pearsons...");
        readPearsons.addActionListener(new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                try {
                    File f = FileDialogUtils.chooseFile("Pearsons file (Yunfan format V2)");
                    if (f != null) {
                        BasicMatrix bm = new DiskResidentBlockMatrix(f.getAbsolutePath());

                        hic.zd.setPearsons(bm);
                    }
                } catch (IOException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }

        });
        extrasMenu.add(readPearsons);

        extrasMenu.add(dumpPearsons);
        if (isInternal) {
            menuBar.add(extrasMenu);
        }
        else {
            resolutionSlider.setMaximum(6);
        }
        return menuBar;
    }


    private JComboBox chrBox1;
    private JComboBox chrBox2;
    private JComboBox displayOptionComboBox;
    private RangeSlider colorRangeSlider;
    private JSlider resolutionSlider;

    //private JScrollPane trackPanelScrollpane;
    TrackPanel trackPanel;
    private HiCRulerPanel rulerPanel2;
    private HeatmapPanel heatmapPanel;
    private HiCRulerPanel rulerPanel1;
    private JCheckBoxMenuItem viewEigenvector;
    private ThumbnailPanel thumbnailPanel;
    private JPanel xPlotPanel;
    private JPanel yPlotPanel;


    enum DisplayOption {
        OBSERVED("Observed"),
        OE("OE"),
        PEARSON("Pearson");
        private String value;

        DisplayOption(String value) {
            this.value = value;
        }

        public String toString() {
            return value;
        }

    }


}
