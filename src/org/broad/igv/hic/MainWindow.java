/*
 * Created by JFormDesigner on Mon Aug 02 22:04:22 EDT 2010
 */

package org.broad.igv.hic;

import java.awt.event.*;
import javax.swing.border.*;
import javax.swing.event.*;

import org.broad.igv.hic.data.*;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import javax.swing.*;

/**
 * @author James Robinson
 */
public class MainWindow extends JFrame {

    int refMaxCount = 500;

    public static int[] zoomBinSizes = {2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 2500, 1000};
    public static final int MAX_ZOOM = 10;

    String genomeId;
    Chromosome chr1;
    Chromosome chr2;
    private int len;
    public Context xContext;
    public Context yContext;
    Dataset dataset;
    MatrixZoomData zd;
    java.util.List<Chromosome> chromosomes;

    public static void main(String[] args) throws IOException {

        final MainWindow mainWindow = new MainWindow();

        mainWindow.pack();
        mainWindow.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        mainWindow.setVisible(true);


    }


    public MainWindow() throws IOException {
        initComponents();


        thumbnailPanel.setMainWindow(this);

        //dataset = AlignmentsParser.read(new File(file));

        updateGenome();


        rangeScale.setValue(refMaxCount);

        heatmapPanel.setSize(500, 500);
        heatmapPanel.setMainWindow(this);

        MouseAdapter mh = new HeatmapMouseHandler();
        heatmapPanel.addMouseListener(mh);
        heatmapPanel.addMouseMotionListener(mh);

        thumbnailPanel.setPreferredSize(new Dimension(100, 100));
        thumbnailPanel.setBinWidth(1);
        thumbnailPanel.setMaxCount(refMaxCount);
        thumbnailPanel.setContext(xContext);

        //2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 2500, 1000
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
            dataset = (new BinDatasetReader(new File(file))).read();
        } else {
            dataset = (new AsciiDatasetReader(genomeId, new File(file))).read();
        }
        setTitle(file);
        refreshChromosomes();
    }

    private void updateGenome() {
        genomeId = (String) genomeBox.getSelectedItem();
        chromosomes = Chromosome.chromosomes.get(genomeId);
        chrBox1.setModel(new DefaultComboBoxModel(chromosomes.toArray()));
        chrBox2.setModel(new DefaultComboBoxModel(chromosomes.toArray()));
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

                chr1 = (Chromosome) chrBox1.getSelectedItem();
                chr2 = (Chromosome) chrBox2.getSelectedItem();
                xContext = new Context(chr1);
                yContext = new Context(chr2);
                rulerPanel2.setFrame(xContext, HiCRulerPanel.Orientation.HORIZONTAL);
                rulerPanel1.setFrame(yContext, HiCRulerPanel.Orientation.VERTICAL);

                int t1 = chr1.getIndex();
                int t2 = chr2.getIndex();

                if (t1 > t2) {
                    Chromosome tmp = chr2;
                    chr2 = chr1;
                    chr1 = tmp;
                }

                Matrix m = dataset.getMatrix(chr1.getIndex(), chr2.getIndex());
                if (m == null) {
                } else {
                    MatrixZoomData thumbnailData = m.getZoomData(0);
                    thumbnailPanel.setZd(thumbnailData);
                    setInitialZoom();
                }


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
        int binWidth = (pixels / nBins);
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
        thumbnailPanel.setZd(dataset.getMatrix(chr1.getIndex(), chr2.getIndex()).getZoomData(thumbnailZoom));
        System.out.println("tz=" + thumbnailZoom);

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

        //int currentBinSize = zoomBinSizes[xContext.getZoom()];
        int newBinSize = zoomBinSizes[newZoom];
        //double ratio = Math.pow(((double) newBinSize) / currentBinSize, 2);
        //int scale = Math.max(1, (int) (ratio * refMaxCount));
        //rangeScale.setValue(scale);

        xContext.setZoom(newZoom, ((double) newBinSize) / heatmapPanel.getBinWidth());
        yContext.setZoom(newZoom, ((double) newBinSize) / heatmapPanel.getBinWidth());

        zd = dataset.getMatrix(chr1.getIndex(), chr2.getIndex()).getZoomData(xContext.getZoom());

        xContext.setVisibleWidth((int) ((double) heatmapPanel.getVisibleRect().width / heatmapPanel.getBinWidth() * zoomBinSizes[xContext.getZoom()]));
        yContext.setVisibleWidth((int) ((double) heatmapPanel.getVisibleRect().height / heatmapPanel.getBinWidth() * zoomBinSizes[xContext.getZoom()]));

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
        xContext.increment(dx);
        yContext.increment(dy);
        repaint();
    }

    private void rangeStateChanged(ChangeEvent e) {
        JSlider slider = (JSlider) e.getSource();
        int maxCount = slider.getValue();

        heatmapPanel.setMaxCount(maxCount);
        thumbnailPanel.setMaxCount(maxCount);
        repaint();
    }


    private void thumbnailPanelMouseClicked(MouseEvent e) {
        int bw = thumbnailPanel.getBinWidth();
        int nBins = getLen() / thumbnailPanel.getZd().getBinSize();
        int effectiveWidth = nBins * bw;

        int newX = (int) (((double) e.getX() / effectiveWidth) * getLen());
        int newY = (int) (((double) e.getY() / effectiveWidth) * getLen());
        center(newX, newY);
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

    private void genomeBoxActionPerformed(ActionEvent e) {
        updateGenome();
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

    class HeatmapMouseHandler extends MouseAdapter {

        boolean isDragging = false;
        private Point lastMousePoint;
        private int lastMousePressedY;
        private int cumulativeDeltaX;
        private int cumulativeDeltaY;


        @Override
        public void mousePressed(final MouseEvent e) {

            isDragging = true;

            //panel = (Container) e.getSource();
            //panel.setCursor(dragCursor);

            lastMousePoint = e.getPoint();
            lastMousePressedY = (int) e.getPoint().getY();
            cumulativeDeltaX = 0;
            cumulativeDeltaY = 0;


        }


        @Override
        public void mouseReleased(final MouseEvent e) {
            if (isDragging) {
                isDragging = false;
            }
            lastMousePoint = null;
            // ((JComponent) e.getSource()).setCursor(getCursor());
        }


        @Override
        final public void mouseDragged(final MouseEvent e) {

            // ((JComponent) e.getSource()).setCursor(IGVMainFrame.handCursor);
            try {

                if (lastMousePoint == null) {
                    lastMousePoint = e.getPoint();
                    return;
                }

                double deltaX = lastMousePoint.getX() - e.getX();
                double deltaY = lastMousePoint.getY() - e.getY();

                // Size i
                int dx = (int) (deltaX * zd.getBinSize() / heatmapPanel.getBinWidth());
                int dy = (int) (deltaY * zd.getBinSize() / heatmapPanel.getBinWidth());

                moveBy(dx, dy);


            } finally {
                lastMousePoint = e.getPoint();    // Always save the last Point
            }
        }

        @Override
        public void mouseClicked(MouseEvent e) {

            if (!e.isPopupTrigger() && (e.getClickCount() > 1)) {
                int currentZoom = xContext.getZoom();
                final int newZoom = e.isAltDown()
                        ? Math.max(currentZoom - 1, 1)
                        : Math.min(11, currentZoom + 1);

                int centerLocationX = xContext.getOrigin() + e.getX() * zd.getBinSize() / heatmapPanel.getBinWidth();
                int centerLocationY = yContext.getOrigin() + e.getY() * zd.getBinSize() / heatmapPanel.getBinWidth();

                setZoom(newZoom, centerLocationX, centerLocationY);

            }
        }

        @Override
        public void mouseMoved(MouseEvent e) {
            if (xContext != null && zd != null) {
                Rectangle visRect = heatmapPanel.getVisibleRect();
                final int binSize = zd.getBinSize();
                int binX = (xContext.getOrigin() / binSize) + (e.getX() - visRect.x) / heatmapPanel.getBinWidth();
                int binY = (yContext.getOrigin() / binSize) + (e.getY() - visRect.y) / heatmapPanel.getBinWidth();
                StringBuffer txt = new StringBuffer();
                txt.append("<html>");
                txt.append(chr1.getName());
                txt.append(":");
                txt.append(String.valueOf((binX - 1) * binSize));
                txt.append("-");
                txt.append(String.valueOf(binX * binSize));
                txt.append("<br>");
                txt.append(chr2.getName());
                txt.append(":");
                txt.append(String.valueOf((binY - 1) * binSize));
                txt.append("-");
                txt.append(String.valueOf(binY * binSize));
                heatmapPanel.setToolTipText(txt.toString());
            }
        }
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        panel2 = new JPanel();
        panel4 = new JPanel();
        genomeBox = new JComboBox();
        chrSelectionPanel = new JPanel();
        chrBox1 = new JComboBox();
        chrBox2 = new JComboBox();
        refreshButton = new JButton();
        zoomComboBox = new JComboBox();
        thumbnailPanel = new ThumbnailPanel();
        panel3 = new JPanel();
        panel5 = new JPanel();
        rulerPanel2 = new HiCRulerPanel();
        panel6 = new JPanel();
        heatmapPanel = new HeatmapPanel();
        rulerPanel1 = new HiCRulerPanel();
        panel1 = new JPanel();
        rangeScale = new JSlider();
        menuBar1 = new JMenuBar();
        fileMenu = new JMenu();
        loadMenuItem = new JMenuItem();

        //======== this ========
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== panel2 ========
        {
            panel2.setLayout(new BorderLayout());

            //======== panel4 ========
            {
                panel4.setLayout(new FlowLayout(FlowLayout.LEFT, 25, 5));

                //---- genomeBox ----
                genomeBox.setModel(new DefaultComboBoxModel(new String[] {
                    "dmel",
                    "hg18"
                }));
                genomeBox.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        genomeBoxActionPerformed(e);
                    }
                });
                panel4.add(genomeBox);

                //======== chrSelectionPanel ========
                {
                    chrSelectionPanel.setBorder(new BevelBorder(BevelBorder.LOWERED));
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

                //---- zoomComboBox ----
                zoomComboBox.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        zoomComboBoxActionPerformed(e);
                    }
                });
                panel4.add(zoomComboBox);

                //---- thumbnailPanel ----
                thumbnailPanel.setMaximumSize(new Dimension(100, 100));
                thumbnailPanel.setMinimumSize(new Dimension(100, 100));
                thumbnailPanel.setPreferredSize(new Dimension(100, 100));
                thumbnailPanel.setBorder(LineBorder.createBlackLineBorder());
                thumbnailPanel.addMouseListener(new MouseAdapter() {
                    @Override
                    public void mouseClicked(MouseEvent e) {
                        thumbnailPanelMouseClicked(e);
                    }
                });
                panel4.add(thumbnailPanel);
            }
            panel2.add(panel4, BorderLayout.NORTH);

            //======== panel3 ========
            {
                panel3.setLayout(new BorderLayout());

                //======== panel5 ========
                {
                    panel5.setLayout(new BorderLayout());

                    //---- rulerPanel2 ----
                    rulerPanel2.setMaximumSize(new Dimension(4000, 50));
                    rulerPanel2.setMinimumSize(new Dimension(1, 50));
                    rulerPanel2.setPreferredSize(new Dimension(1, 50));
                    rulerPanel2.setBorder(LineBorder.createBlackLineBorder());
                    panel5.add(rulerPanel2, BorderLayout.CENTER);

                    //======== panel6 ========
                    {
                        panel6.setMaximumSize(new Dimension(50, 50));
                        panel6.setMinimumSize(new Dimension(50, 50));
                        panel6.setPreferredSize(new Dimension(50, 50));
                        panel6.setLayout(null);
                    }
                    panel5.add(panel6, BorderLayout.WEST);
                }
                panel3.add(panel5, BorderLayout.NORTH);

                //---- heatmapPanel ----
                heatmapPanel.setBorder(LineBorder.createBlackLineBorder());
                heatmapPanel.setMaximumSize(new Dimension(800, 800));
                heatmapPanel.setMinimumSize(new Dimension(800, 800));
                heatmapPanel.setPreferredSize(new Dimension(800, 800));
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
                rulerPanel1.setBorder(LineBorder.createBlackLineBorder());
                rulerPanel1.setMinimumSize(new Dimension(50, 1));
                panel3.add(rulerPanel1, BorderLayout.WEST);
            }
            panel2.add(panel3, BorderLayout.CENTER);

            //======== panel1 ========
            {
                panel1.setLayout(new FlowLayout());

                //---- rangeScale ----
                rangeScale.setMajorTickSpacing(100);
                rangeScale.setPaintTicks(true);
                rangeScale.setPaintLabels(true);
                rangeScale.setPreferredSize(new Dimension(300, 52));
                rangeScale.setValue(500);
                rangeScale.setMaximum(1000);
                rangeScale.addChangeListener(new ChangeListener() {
                    public void stateChanged(ChangeEvent e) {
                        rangeStateChanged(e);
                    }
                });
                panel1.add(rangeScale);
            }
            panel2.add(panel1, BorderLayout.SOUTH);
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
            }
            menuBar1.add(fileMenu);
        }
        contentPane.add(menuBar1, BorderLayout.NORTH);
        pack();
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel panel2;
    private JPanel panel4;
    private JComboBox genomeBox;
    private JPanel chrSelectionPanel;
    private JComboBox chrBox1;
    private JComboBox chrBox2;
    private JButton refreshButton;
    private JComboBox zoomComboBox;
    private ThumbnailPanel thumbnailPanel;
    private JPanel panel3;
    private JPanel panel5;
    private HiCRulerPanel rulerPanel2;
    private JPanel panel6;
    private HeatmapPanel heatmapPanel;
    private HiCRulerPanel rulerPanel1;
    private JPanel panel1;
    private JSlider rangeScale;
    private JMenuBar menuBar1;
    private JMenu fileMenu;
    private JMenuItem loadMenuItem;
    // JFormDesigner - End of variables declaration  //GEN-END:variables


}
