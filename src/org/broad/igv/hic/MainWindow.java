/*
 * Created by JFormDesigner on Mon Aug 02 22:04:22 EDT 2010
 */

package org.broad.igv.hic;

import java.awt.event.*;
import javax.swing.border.*;
import javax.swing.event.*;

import org.broad.igv.hic.data.*;
import sun.security.krb5.internal.crypto.Crc32CksumType;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.util.*;
import javax.swing.*;

/**
 * @author James Robinson
 */
public class MainWindow extends JFrame {

    //String file = "test/data/test.summary.binned.sorted.hic";
    //String file = "data/GSE18199_RAW/GSM455140_428EGAAXX.8.maq.hic.summary.binned.txt";
    //String file = "data/GSE18199_RAW/GSM455140_428EGAAXX.8.maq.hic.summary.binned.compressed.hic";

    static int refMaxCount = 30;
    public static int[] zoomBinSizes = {2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 2500, 1000};
    public static final int MAX_ZOOM = 10;

    String genomeId;
    Chromosome chr1;
    Chromosome chr2;
    int len;
    public Context context;
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

        context = new Context();

        rangeScale.setValue(refMaxCount);

        heatmapPanel.setSize(500, 500);
        heatmapPanel.setMainWindow(this);

        MouseAdapter mh = new HeatmapMouseHandler();
        heatmapPanel.addMouseListener(mh);
        heatmapPanel.addMouseMotionListener(mh);

        thumbnailPanel.setPreferredSize(new Dimension(100, 100));
        thumbnailPanel.setBinWidth(1);
        thumbnailPanel.setMaxCount((int) (6.25 * refMaxCount));
        thumbnailPanel.setContext(context);

        Hashtable<Integer, JLabel> zoomLabels = new Hashtable();
        zoomLabels.put(1, new JLabel("1 MB"));
        zoomLabels.put(4, new JLabel("100 KB"));
        zoomLabels.put(7, new JLabel("10 KB"));
        zoomLabels.put(10, new JLabel("1 KB"));
        zoomSlider.setLabelTable(zoomLabels);
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

        context.setLenX(chr1.getSize());
        context.setLenY(chr2.getSize());
        len = Math.max(context.getLenX(), context.getLenY());
        int pixels = heatmapPanel.getVisibleRect().width;
        int maxNBins = pixels / 2;

        // Main panel
        int bp_bin = len / maxNBins;
        int initialZoom = zoomBinSizes.length - 1;
        for (int z = 1; z < zoomBinSizes.length; z++) {
            if (zoomBinSizes[z] < bp_bin) {
                initialZoom = z - 1;
                break;
            }
        }
        int nBins = len / zoomBinSizes[initialZoom] + 1;
        int binWidth = (pixels / nBins);
        heatmapPanel.setBinWidth(binWidth);


        // Thumbnail
        pixels = thumbnailPanel.getVisibleRect().width;
        maxNBins = pixels;
        bp_bin = len / maxNBins;
        int thumbnailZoom = zoomBinSizes.length - 1;
        for (int z = 1; z < zoomBinSizes.length; z++) {
            if (zoomBinSizes[z] < bp_bin) {
                thumbnailZoom = z - 1;
                break;
            }
        }
        nBins = len / zoomBinSizes[thumbnailZoom] + 1;
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
        int centerLocationX = context.getOriginX() + (int) (visibleRect.getWidth() * binSize / (2 * heatmapPanel.getBinWidth()));
        int centerLocationY = context.getOriginY() + (int) (visibleRect.getWidth() * binSize / (2 * heatmapPanel.getBinWidth()));
        setZoom(zoom, centerLocationX, centerLocationY);
    }

    public void setZoom(int newZoom, int centerLocationX, int centerLocationY) {

        if (newZoom < 1 || newZoom > MAX_ZOOM) return;

        int newBinSize = zoomBinSizes[newZoom];
        double ratio = ((double) newBinSize) / zoomBinSizes[1];
        int scale = Math.max(1, (int) (ratio * refMaxCount));
        rangeScale.setValue(scale);

        context.setZoom(newZoom);

        zd = dataset.getMatrix(chr1.getIndex(), chr2.getIndex()).getZoomData(context.getZoom());

        context.setVisibleWidth((int) ((double) heatmapPanel.getVisibleRect().width / heatmapPanel.getBinWidth() * zoomBinSizes[context.getZoom()]));
        context.setVisibleHeight((int) ((double) heatmapPanel.getVisibleRect().height / heatmapPanel.getBinWidth() * zoomBinSizes[context.getZoom()]));

        zoomSlider.setValue(newZoom);

        center(centerLocationX, centerLocationY);

    }

    public void center(int centerLocationX, int centerLocationY) {

        int binSize = zd.getBinSize();
        int w = (int) (heatmapPanel.getVisibleRect().getWidth() * binSize / heatmapPanel.getBinWidth());
        context.setOrigin(((int) (centerLocationX - w / 2)), (int) (centerLocationY - w / 2));
        heatmapPanel.repaint();
        thumbnailPanel.repaint();
    }


    public void moveBy(int dx, int dy) {
        context.increment(dx, dy);
        heatmapPanel.repaint();
        thumbnailPanel.repaint();
    }

    private void rangeStateChanged(ChangeEvent e) {
        JSlider slider = (JSlider) e.getSource();
        int maxCount = Math.max(1, slider.getValue());

        heatmapPanel.setMaxCount(maxCount); //getZoomMaxCount());
        heatmapPanel.repaint();
    }

    private void zoomSliderStateChanged(ChangeEvent e) {

        int newZoom = zoomSlider.getValue();
        if (newZoom == context.getZoom()) return;
        setZoom(newZoom);
    }

    public int getZoomMaxCount() {

        int factor = (int) Math.pow((double) zoomBinSizes[context.getZoom()] / zoomBinSizes[1], 2);
        return (int) Math.max(1, (refMaxCount * factor));
    }

    private void thumbnailPanelMouseClicked(MouseEvent e) {
        int bw = thumbnailPanel.getBinWidth();
        int nBins = len / thumbnailPanel.getZd().getBinSize();
        int effectiveWidth = nBins * bw;

        int newX = (int) (((double) e.getX() / effectiveWidth) * len);
        int newY = (int) (((double) e.getY() / effectiveWidth) * len);
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
                int currentZoom = context.getZoom();
                final int newZoom = e.isAltDown()
                        ? Math.max(currentZoom - 1, 1)
                        : Math.min(11, currentZoom + 1);

                int centerLocationX = context.getOriginX() + e.getX() * zd.getBinSize() / heatmapPanel.getBinWidth();
                int centerLocationY = context.getOriginY() + e.getY() * zd.getBinSize() / heatmapPanel.getBinWidth();

                setZoom(newZoom, centerLocationX, centerLocationY);

            }
        }

        @Override
        public void mouseMoved(MouseEvent e) {
            if (context != null && zd != null) {
                Rectangle visRect = heatmapPanel.getVisibleRect();
                final int binSize = zd.getBinSize();
                int binX = (context.getOriginX() / binSize) + (e.getX() - visRect.x) / heatmapPanel.getBinWidth();
                int binY = (context.getOriginY() / binSize) + (e.getY() - visRect.y) / heatmapPanel.getBinWidth();
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
        zoomSlider = new JSlider();
        heatmapPanel = new HeatmapPanel();
        panel1 = new JPanel();
        rangeScale = new JSlider();
        thumbnailPanel = new ThumbnailPanel();
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

                //---- zoomSlider ----
                zoomSlider.setMinimum(1);
                zoomSlider.setMaximum(10);
                zoomSlider.setMajorTickSpacing(3);
                zoomSlider.setPaintLabels(true);
                zoomSlider.setValue(1);
                zoomSlider.setSnapToTicks(true);
                zoomSlider.setName("Zoom slider");
                zoomSlider.setPaintTicks(true);
                zoomSlider.setMinorTickSpacing(1);
                zoomSlider.addChangeListener(new ChangeListener() {
                    public void stateChanged(ChangeEvent e) {
                        zoomSliderStateChanged(e);
                    }
                });
                panel4.add(zoomSlider);
            }
            panel2.add(panel4, BorderLayout.NORTH);

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
            panel2.add(heatmapPanel, BorderLayout.CENTER);

            //======== panel1 ========
            {
                panel1.setLayout(new FlowLayout());

                //---- rangeScale ----
                rangeScale.setMaximum(500);
                rangeScale.setMinorTickSpacing(10);
                rangeScale.setMajorTickSpacing(50);
                rangeScale.setPaintTicks(true);
                rangeScale.setPaintLabels(true);
                rangeScale.setPreferredSize(new Dimension(300, 52));
                rangeScale.addChangeListener(new ChangeListener() {
                    public void stateChanged(ChangeEvent e) {
                        rangeStateChanged(e);
                    }
                });
                panel1.add(rangeScale);

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
                panel1.add(thumbnailPanel);
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
    private JSlider zoomSlider;
    private HeatmapPanel heatmapPanel;
    private JPanel panel1;
    private JSlider rangeScale;
    private ThumbnailPanel thumbnailPanel;
    private JMenuBar menuBar1;
    private JMenu fileMenu;
    private JMenuItem loadMenuItem;
    // JFormDesigner - End of variables declaration  //GEN-END:variables


}
