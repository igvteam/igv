/*
 * Created by JFormDesigner on Mon Aug 02 22:04:22 EDT 2010
 */

package org.broad.igv.hic;

import java.awt.event.*;
import javax.swing.border.*;
import javax.swing.event.*;

import com.jidesoft.swing.*;

import org.broad.igv.hic.data.*;
import org.broad.igv.ui.util.IconFactory;
import org.broad.tribble.util.SeekableStream;
import org.broad.tribble.util.SeekableStreamFactory;

import java.awt.*;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import javax.swing.*;

/**
 * @author James Robinson
 */
public class MainWindow extends JFrame {

    public static final int MIN_BIN_WIDTH = 2;
    public static Cursor fistCursor;
    //int refMaxCount = 500;

    public static int[] zoomBinSizes = {2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 2500, 1000};
    public static final int MAX_ZOOM = 10;

    private int len;
    public Context xContext;
    public Context yContext;
    Dataset dataset;
    MatrixZoomData zd;
    Chromosome[] chromosomes;

    public static void main(String[] args) throws IOException {

        final MainWindow mainWindow = new MainWindow();
        mainWindow.setVisible(true);
        mainWindow.setSize(870, 870);


    }


    public void createCursors() {
        BufferedImage handImage = new BufferedImage(32, 32, BufferedImage.TYPE_INT_ARGB);

        // Make backgroun transparent
        Graphics2D g = handImage.createGraphics();
        g.setComposite(AlphaComposite.getInstance(
                AlphaComposite.CLEAR, 0.0f));
        Rectangle2D.Double rect = new Rectangle2D.Double(0, 0, 32, 32);
        g.fill(rect);

        // Draw hand image in middle
        g = handImage.createGraphics();
        g.drawImage(IconFactory.getInstance().getIcon(IconFactory.IconID.FIST).getImage(), 0, 0, null);
        MainWindow.fistCursor = getToolkit().createCustomCursor(handImage, new Point(8, 6), "Move");
    }


    public MainWindow() throws IOException {

        initComponents();

        createCursors();

        thumbnailPanel.setMainWindow(this);

        int initialMaxCount = 50000;
        rangeScale.setValue(initialMaxCount);

        heatmapPanel.setSize(500, 500);
        heatmapPanel.setMainWindow(this);

        thumbnailPanel.setPreferredSize(new Dimension(100, 100));
        thumbnailPanel.setBinWidth(1);

        //2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 2500, 1000
        // TODO -- these should be read from the data file  (zd.binSize)
        ZoomLabel[] zooms = new ZoomLabel[]{
                new ZoomLabel("2.5 mb", 0),
                new ZoomLabel("1   mb", 1),
                new ZoomLabel("500 kb", 2),
                new ZoomLabel("250 kb", 3),
                new ZoomLabel("100 kb", 4),
                new ZoomLabel("50  kb", 5),
                new ZoomLabel("25  kb", 6),
                new ZoomLabel("10  kb", 7),
                new ZoomLabel("5   kb", 8),
                new ZoomLabel("2.5 kb", 9),
                new ZoomLabel("1   kb", 10)};
        zoomComboBox.setModel(new DefaultComboBoxModel(zooms));

        pack();
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);


    }

    class ZoomLabel {
        String label;
        int zoom;

        ZoomLabel(String label, int zoom) {
            this.label = label;
            this.zoom = zoom;
        }

        public String toString() {
            return label;
        }
    }

    private void load(String file) throws IOException {
        if (file.endsWith("hic")) {
            SeekableStream ss = SeekableStreamFactory.getStreamFor(file);
            dataset = (new DatasetReader(ss)).read();
            chromosomes = dataset.getChromosomes();
            chrBox1.setModel(new DefaultComboBoxModel(chromosomes));
            chrBox2.setModel(new DefaultComboBoxModel(chromosomes));
        } else {
            // error -- unknown file type
        }
        setTitle(file);
        refreshChromosomes();
    }

    private void refreshChromosomes() {

        SwingWorker worker = new SwingWorker() {
            @Override
            protected void done() {
                getGlassPane().setVisible(false);
                getContentPane().repaint();
            }

            @Override
            protected Object doInBackground() throws Exception {
                getGlassPane().setVisible(true);
                getGlassPane().setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));

                Chromosome chr1 = (Chromosome) chrBox1.getSelectedItem();
                Chromosome chr2 = (Chromosome) chrBox2.getSelectedItem();

                int t1 = chr1.getIndex();
                int t2 = chr2.getIndex();

                if (t1 > t2) {
                    Chromosome tmp = chr2;
                    chr2 = chr1;
                    chr1 = tmp;
                }

                xContext = new Context(chr2);
                yContext = new Context(chr1);
                rulerPanel2.setFrame(xContext, HiCRulerPanel.Orientation.HORIZONTAL);
                rulerPanel1.setFrame(yContext, HiCRulerPanel.Orientation.VERTICAL);

                Matrix m = dataset.getMatrix(chr1, chr2);
                if (m == null) {
                } else {
                    setInitialZoom();
                }

                heatmapPanel.clearTileCache();


                Image thumbnail = heatmapPanel.getThumbnailImage(zd, thumbnailPanel.getWidth(), thumbnailPanel.getHeight());
                thumbnailPanel.setImage(thumbnail);

                return null;
            }

        };

        worker.execute();

    }

    private void setInitialZoom() {

        setLen(Math.max(xContext.getChrLength(), yContext.getChrLength()));
        int pixels = heatmapPanel.getVisibleRect().width;
        int maxNBins = pixels / 2;

        // Main panel
        int bp_bin = getLen() / maxNBins;
        int initialZoom = zoomBinSizes.length - 1;
        for (int z = 1; z < zoomBinSizes.length; z++) {
            if (zoomBinSizes[z] < bp_bin) {
                initialZoom = z - 1;
                break;
            }
        }
        int nBins = getLen() / zoomBinSizes[initialZoom] + 1;
        int binWidth = Math.max(MIN_BIN_WIDTH, pixels / nBins);
        heatmapPanel.setBinWidth(binWidth);


        // Thumbnail
        pixels = thumbnailPanel.getVisibleRect().width;
        maxNBins = pixels;
        bp_bin = getLen() / maxNBins;
        int thumbnailZoom = zoomBinSizes.length - 1;
        for (int z = 1; z < zoomBinSizes.length; z++) {
            if (zoomBinSizes[z] < bp_bin) {
                thumbnailZoom = z - 1;
                break;
            }
        }
        nBins = getLen() / zoomBinSizes[thumbnailZoom] + 1;
        binWidth = (pixels / nBins);
        thumbnailPanel.setBinWidth(binWidth);
        setZoom(initialZoom);

    }

    public void setZoom(int zoom) {

        if (zoom < 1 || zoom > MAX_ZOOM) return;

        Rectangle visibleRect = heatmapPanel.getVisibleRect();
        int binSize = zoomBinSizes[zoom];
        int centerLocationX = xContext.getOrigin() + (int) (visibleRect.getWidth() * binSize / (2 * heatmapPanel.getBinWidth()));
        int centerLocationY = yContext.getOrigin() + (int) (visibleRect.getWidth() * binSize / (2 * heatmapPanel.getBinWidth()));
        setZoom(zoom, centerLocationX, centerLocationY);
    }

    public void setZoom(int newZoom, int centerLocationX, int centerLocationY) {

        if (newZoom < 1 || newZoom > MAX_ZOOM) return;

        Chromosome chr1 = xContext.getChromosome();
        Chromosome chr2 = yContext.getChromosome();
        zd = dataset.getMatrix(chr1, chr2).getZoomData(newZoom);

        int newBinSize = zd.getBinSize();

        // Adjust bin width to fill space, if possible.  Bins are square at all times
        int nBinsX = chr1.getSize() / newBinSize + 1;
        int binWidthX = Math.max(MIN_BIN_WIDTH, heatmapPanel.getWidth() / nBinsX);
        int nBinsY = chr2.getSize() / newBinSize + 1;
        int binWidthY = Math.max(MIN_BIN_WIDTH, heatmapPanel.getHeight() / nBinsY);
        int newBinWidth = Math.min(binWidthX, binWidthY);
        heatmapPanel.setBinWidth(newBinWidth);

        // Scale in basepairs per screen pixel
        double scale = (double) newBinSize / newBinWidth;
        xContext.setZoom(newZoom, scale);
        yContext.setZoom(newZoom, scale);


        xContext.setVisibleWidth((int) (scale * heatmapPanel.getWidth()));
        yContext.setVisibleWidth((int) (scale * heatmapPanel.getHeight()));

        zoomComboBox.setSelectedIndex(newZoom);

        center(centerLocationX, centerLocationY);

        repaint();

    }

    public void center(int centerLocationX, int centerLocationY) {

        int binSize = zd.getBinSize();
        int w = (int) (heatmapPanel.getVisibleRect().getWidth() * binSize / heatmapPanel.getBinWidth());
        xContext.setOrigin((int) (centerLocationX - w / 2));
        yContext.setOrigin((int) (centerLocationY - w / 2));
        //heatmapPanel.repaint();
        //thumbnailPanel.repaint();
        repaint();
    }


    public void moveBy(int dx, int dy) {
        xContext.moveBy(dx);
        yContext.moveBy(dy);
        repaint();
    }

    private void rangeStateChanged(ChangeEvent e) {
        JSlider slider = (JSlider) e.getSource();
        int maxCount = slider.getValue();

        heatmapPanel.setMaxCount(maxCount);
        repaint();
    }


    private void heatmapPanelMouseDragged(MouseEvent e) {
        // TODO add your code here
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


    public int getLen() {
        return len;
    }

    public void setLen(int len) {
        this.len = len;
    }

    private void zoomComboBoxActionPerformed(ActionEvent e) {
        Object selected = zoomComboBox.getSelectedItem();
        if (selected != null) {
            int newZoom = ((ZoomLabel) selected).zoom;
            if (newZoom == xContext.getZoom()) return;
            setZoom(newZoom);
            repaint();
        }

    }

    private void exitActionPerformed(ActionEvent e) {
        setVisible(false);
        dispose();
        System.exit(0);
    }

    private void loadDmelDatasetActionPerformed(ActionEvent e) {
        try {
            load("http://iwww.broadinstitute.org/igvdata/hic/selected_formatted.hic");
        } catch (IOException e1) {
            JOptionPane.showMessageDialog(this, "Error loading data: " + e1.getMessage());
        }
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        panel2 = new JPanel();
        panel4 = new JPanel();
        chrSelectionPanel = new JPanel();
        chrBox1 = new JComboBox();
        chrBox2 = new JComboBox();
        refreshButton = new JButton();
        panel7 = new JPanel();
        label2 = new JLabel();
        zoomComboBox = new JComboBox();
        panel1 = new JPanel();
        panel9 = new JPanel();
        label1 = new JLabel();
        rangeScale = new JSlider();
        panel3 = new JPanel();
        panel5 = new JPanel();
        spacerLeft = new JPanel();
        rulerPanel2 = new HiCRulerPanel();
        spacerRight = new JPanel();
        heatmapPanel = new HeatmapPanel();
        rulerPanel1 = new HiCRulerPanel();
        panel8 = new JPanel();
        thumbnailPanel = new ThumbnailPanel();
        menuBar1 = new JMenuBar();
        fileMenu = new JMenu();
        loadMenuItem = new JMenuItem();
        loadFromURL = new JMenuItem();
        loadDmelDataset = new JMenuItem();
        exit = new JMenuItem();

        //======== this ========
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== panel2 ========
        {
            panel2.setLayout(new BorderLayout());

            //======== panel4 ========
            {
                panel4.setBorder(new MatteBorder(1, 1, 1, 1, Color.black));
                panel4.setLayout(new FlowLayout(FlowLayout.LEFT, 25, 5));

                //======== chrSelectionPanel ========
                {
                    chrSelectionPanel.setBorder(new MatteBorder(1, 1, 1, 1, Color.black));
                    chrSelectionPanel.setLayout(new FlowLayout());
                    chrSelectionPanel.add(chrBox1);
                    chrSelectionPanel.add(chrBox2);

                    //---- refreshButton ----
                    refreshButton.setText("Refresh");
                    refreshButton.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent e) {
                            refreshButtonActionPerformed(e);
                        }
                    });
                    chrSelectionPanel.add(refreshButton);
                }
                panel4.add(chrSelectionPanel);

                //======== panel7 ========
                {
                    panel7.setBorder(new MatteBorder(1, 1, 1, 1, Color.black));
                    panel7.setLayout(new FlowLayout());

                    //---- label2 ----
                    label2.setText("Resolution");
                    panel7.add(label2);

                    //---- zoomComboBox ----
                    zoomComboBox.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent e) {
                            zoomComboBoxActionPerformed(e);
                        }
                    });
                    panel7.add(zoomComboBox);
                }
                panel4.add(panel7);

                //======== panel1 ========
                {
                    panel1.setLayout(new FlowLayout());

                    //======== panel9 ========
                    {
                        panel9.setBorder(null);
                        panel9.setLayout(new FlowLayout());

                        //---- label1 ----
                        label1.setText("Scale");
                        panel9.add(label1);

                        //---- rangeScale ----
                        rangeScale.setMajorTickSpacing(50000);
                        rangeScale.setPaintTicks(true);
                        rangeScale.setPaintLabels(true);
                        rangeScale.setPreferredSize(new Dimension(150, 52));
                        rangeScale.setValue(500);
                        rangeScale.setMaximum(100000);
                        rangeScale.setMinorTickSpacing(10000);
                        rangeScale.addChangeListener(new ChangeListener() {
                            public void stateChanged(ChangeEvent e) {
                                rangeStateChanged(e);
                            }
                        });
                        panel9.add(rangeScale);
                    }
                    panel1.add(panel9);
                }
                panel4.add(panel1);
            }
            panel2.add(panel4, BorderLayout.NORTH);

            //======== panel3 ========
            {
                panel3.setLayout(new BorderLayout());

                //======== panel5 ========
                {
                    panel5.setLayout(new BorderLayout());

                    //======== spacerLeft ========
                    {
                        spacerLeft.setMaximumSize(new Dimension(50, 50));
                        spacerLeft.setMinimumSize(new Dimension(50, 50));
                        spacerLeft.setPreferredSize(new Dimension(50, 50));
                        spacerLeft.setLayout(null);
                    }
                    panel5.add(spacerLeft, BorderLayout.WEST);

                    //---- rulerPanel2 ----
                    rulerPanel2.setMaximumSize(new Dimension(4000, 50));
                    rulerPanel2.setMinimumSize(new Dimension(1, 50));
                    rulerPanel2.setPreferredSize(new Dimension(1, 50));
                    rulerPanel2.setBorder(null);
                    panel5.add(rulerPanel2, BorderLayout.CENTER);

                    //======== spacerRight ========
                    {
                        spacerRight.setMinimumSize(new Dimension(120, 0));
                        spacerRight.setPreferredSize(new Dimension(120, 0));
                        spacerRight.setMaximumSize(new Dimension(120, 32767));
                        spacerRight.setLayout(null);
                    }
                    panel5.add(spacerRight, BorderLayout.EAST);
                }
                panel3.add(panel5, BorderLayout.NORTH);

                //---- heatmapPanel ----
                heatmapPanel.setBorder(LineBorder.createBlackLineBorder());
                heatmapPanel.setMaximumSize(new Dimension(500, 500));
                heatmapPanel.setMinimumSize(new Dimension(500, 500));
                heatmapPanel.setPreferredSize(new Dimension(500, 500));
                heatmapPanel.addMouseMotionListener(new MouseMotionAdapter() {
                    @Override
                    public void mouseDragged(MouseEvent e) {
                        heatmapPanelMouseDragged(e);
                    }
                });
                panel3.add(heatmapPanel, BorderLayout.CENTER);

                //---- rulerPanel1 ----
                rulerPanel1.setMaximumSize(new Dimension(50, 4000));
                rulerPanel1.setPreferredSize(new Dimension(50, 500));
                rulerPanel1.setBorder(null);
                rulerPanel1.setMinimumSize(new Dimension(50, 1));
                panel3.add(rulerPanel1, BorderLayout.WEST);

                //======== panel8 ========
                {
                    panel8.setMaximumSize(new Dimension(120, 100));
                    panel8.setBorder(new EmptyBorder(0, 10, 0, 0));
                    panel8.setLayout(new FlowLayout());

                    //---- thumbnailPanel ----
                    thumbnailPanel.setMaximumSize(new Dimension(100, 100));
                    thumbnailPanel.setMinimumSize(new Dimension(100, 100));
                    thumbnailPanel.setPreferredSize(new Dimension(100, 100));
                    thumbnailPanel.setBorder(LineBorder.createBlackLineBorder());
                    panel8.add(thumbnailPanel);
                }
                panel3.add(panel8, BorderLayout.EAST);
            }
            panel2.add(panel3, BorderLayout.CENTER);
        }
        contentPane.add(panel2, BorderLayout.CENTER);

        //======== menuBar1 ========
        {

            //======== fileMenu ========
            {
                fileMenu.setText("File");

                //---- loadMenuItem ----
                loadMenuItem.setText("Load...");
                loadMenuItem.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        loadMenuItemActionPerformed(e);
                    }
                });
                fileMenu.add(loadMenuItem);

                //---- loadFromURL ----
                loadFromURL.setText("Load from URL ...");
                loadFromURL.setName("loadFromURL");
                loadFromURL.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        loadFromURLActionPerformed(e);
                    }
                });
                fileMenu.add(loadFromURL);
                fileMenu.addSeparator();

                //---- loadDmelDataset ----
                loadDmelDataset.setText("https://iwww.broadinstitute.org/igvdata/hic/selected_formatted.hic");
                loadDmelDataset.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        loadDmelDatasetActionPerformed(e);
                    }
                });
                fileMenu.add(loadDmelDataset);
                fileMenu.addSeparator();

                //---- exit ----
                exit.setText("Exit");
                exit.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        exitActionPerformed(e);
                    }
                });
                fileMenu.add(exit);
            }
            menuBar1.add(fileMenu);
        }
        contentPane.add(menuBar1, BorderLayout.NORTH);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel panel2;
    private JPanel panel4;
    private JPanel chrSelectionPanel;
    private JComboBox chrBox1;
    private JComboBox chrBox2;
    private JButton refreshButton;
    private JPanel panel7;
    private JLabel label2;
    private JComboBox zoomComboBox;
    private JPanel panel1;
    private JPanel panel9;
    private JLabel label1;
    private JSlider rangeScale;
    private JPanel panel3;
    private JPanel panel5;
    private JPanel spacerLeft;
    private HiCRulerPanel rulerPanel2;
    private JPanel spacerRight;
    private HeatmapPanel heatmapPanel;
    private HiCRulerPanel rulerPanel1;
    private JPanel panel8;
    private ThumbnailPanel thumbnailPanel;
    private JMenuBar menuBar1;
    private JMenu fileMenu;
    private JMenuItem loadMenuItem;
    private JMenuItem loadFromURL;
    private JMenuItem loadDmelDataset;
    private JMenuItem exit;
    // JFormDesigner - End of variables declaration  //GEN-END:variables


}
