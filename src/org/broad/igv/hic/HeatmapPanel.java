package org.broad.igv.hic;

import com.jidesoft.swing.JidePopupMenu;
import org.broad.igv.Globals;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.hic.data.MatrixZoomData;
import org.broad.igv.hic.track.HiCFragmentAxis;
import org.broad.igv.hic.track.HiCGridAxis;
import org.broad.igv.util.ObjectCache;

import javax.swing.*;
import javax.swing.text.NumberFormatter;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.io.*;
import java.text.NumberFormat;

/**
 * @author jrobinso
 * @date Aug 2, 2010
 */
public class HeatmapPanel extends JComponent implements Serializable {


    NumberFormat formatter = NumberFormat.getInstance();

    enum DragMode {NONE, PAN, ZOOM}

    ;

    MainWindow mainWindow;
    private HiC hic;
    private int lastCursorX;

    /**
     * Image tile width in pixels
     */
    private int imageTileWidth = 500;

    ObjectCache<String, ImageTile> tileCache = new ObjectCache<String, ImageTile>(10);
    private Rectangle zoomRectangle;

    /**
     * Chromosome boundaries in kbases for whole genome view.
     */
    private int[] chromosomeBoundaries;

    HeatmapRenderer renderer;

    public HeatmapPanel(MainWindow mainWindow, HiC hic) {
        this.mainWindow = mainWindow;
        this.hic = hic;
        renderer = new HeatmapRenderer(mainWindow, hic);
        final HeatmapMouseHandler mouseHandler = new HeatmapMouseHandler();
        addMouseListener(mouseHandler);
        addMouseMotionListener(mouseHandler);
    }


    public void setObservedRange(int min, int max) {
        renderer.setObservedRange(min, max);
        clearTileCache();
        repaint();

    }

    public void setChromosomeBoundaries(int[] chromosomeBoundaries) {
        this.chromosomeBoundaries = chromosomeBoundaries;
    }

    @Override
    protected void paintComponent(Graphics g) {

        ((Graphics2D) g).setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        Rectangle clipBounds = g.getClipBounds();
        g.clearRect(clipBounds.x, clipBounds.y, clipBounds.width, clipBounds.height);

        if (hic != null && hic.zd != null) {
            if (hic.xContext == null)
                return;

            // Same scale used for X & Y (square pixels)
            double scalefactor = hic.xContext.getScaleFactor();

            // Bounds in bins
            int binLeft = hic.xContext.getBinOrigin();
            int bRight = binLeft + (int) (getBounds().width / scalefactor);
            int bTop = hic.yContext.getBinOrigin();
            int bBottom = bTop + (int) (getBounds().height / scalefactor);

            // tile coordinates
            int tLeft = binLeft / imageTileWidth;
            int tRight = bRight / imageTileWidth;
            int tTop = bTop / imageTileWidth;
            int tBottom = bBottom / imageTileWidth;

            // Pixels per bin -- used to scale image


            for (int tileRow = tTop; tileRow <= tBottom; tileRow++) {
                for (int tileColumn = tLeft; tileColumn <= tRight; tileColumn++) {
                    ImageTile tile = getImageTile(tileRow, tileColumn, hic.getDisplayOption());
                    if (tile != null) {
                        if (scalefactor > 0.99 && scalefactor < 1.01) scalefactor = 1;


                        int pxSrcOffset = tile.bLeft - binLeft;
                        int pxDestOffset = (int) (pxSrcOffset * scalefactor);
                        int pySrcOffset = tile.bTop - bTop;
                        int pyDestOffset = (int) (pySrcOffset * scalefactor);

                        int imageWidth = tile.image.getWidth(null);
                        int imageHeight = tile.image.getHeight(null);
                        int xSrc1, xSrc2, ySrc1, ySrc2, xDest1, yDest1, xDest2, yDest2;
                        if (pxSrcOffset < 0) {
                            xSrc1 = -pxSrcOffset;
                            xSrc2 = imageWidth;
                            xDest1 = 0;
                        } else {
                            xSrc1 = 0;
                            int delta = getWidth() - pxSrcOffset;
                            xSrc2 = Math.min(imageWidth, delta);
                            xDest1 = pxDestOffset;
                        }
                        if (pySrcOffset < 0) {
                            ySrc1 = -pySrcOffset;
                            ySrc2 = imageHeight;
                            yDest1 = 0;
                        } else {
                            ySrc1 = 0;
                            int delta = getHeight() - pySrcOffset;
                            ySrc2 = Math.min(imageHeight, delta);
                            yDest1 = pyDestOffset;
                        }

                        xDest2 = (int) (xDest1 + (xSrc2 - xSrc1) * scalefactor);
                        yDest2 = (int) (yDest1 + (ySrc2 - ySrc1) * scalefactor);
                        if (xDest2 <= xDest1 || yDest2 <= yDest1) continue;

                        g.drawImage(tile.image,
                                xDest1, yDest1, xDest2, yDest2,
                                xSrc1, ySrc1, xSrc2, ySrc2,
                                null);
                    }

                }
            }

            boolean isWholeGenome = (hic.xContext.getChromosome().getName().equals("All") &&
                    hic.yContext.getChromosome().getName().equals("All"));


            // Draw grid
            if (isWholeGenome) {
                Color color = g.getColor();
                g.setColor(Color.lightGray);

                Chromosome[] chromosomes = hic.getChromosomes();
                // Index 0 is whole genome
                int xGenomeCoord = 0;
                for (int i = 1; i < chromosomes.length; i++) {
                    Chromosome c = chromosomes[i];
                    xGenomeCoord += (c.getLength() / 1000);
                    int xBin = hic.zd.getxGridAxis().getBinNumberForGenomicPosition(xGenomeCoord);
                    int x = (int) (xBin * hic.xContext.getScaleFactor());
                    g.drawLine(x, 0, x, getHeight());
                }
                int yGenomeCoord = 0;
                for (int i = 1; i < chromosomes.length; i++) {
                    Chromosome c = chromosomes[i];
                    yGenomeCoord += (c.getLength() / 1000);
                    int yBin = hic.zd.getyGridAxis().getBinNumberForGenomicPosition(yGenomeCoord);
                    int y = (int) (yBin * hic.yContext.getScaleFactor());
                    g.drawLine(0, y, getWidth(), y);
                }

                g.setColor(color);
            }

            int cursorX = hic.getCursorX();
            if (cursorX > 0) {
                g.setColor(MainWindow.RULER_LINE_COLOR);
                g.drawLine(cursorX, 0, cursorX, getHeight());
            }


        }

        if (zoomRectangle != null) {
            ((Graphics2D) g).draw(zoomRectangle);
        }
    }

    public Image getThumbnailImage(MatrixZoomData zd0, int tw, int th, MainWindow.DisplayOption displayOption) {

//        int maxBinCountX = (hic.xContext.getChrLength() - hic.xContext.getOrigin()) / zd.getBinSize() + 1;
//        int maxBinCountY = (hic.yContext.getChrLength() - hic.yContext.getOrigin()) / zd.getBinSize() + 1;
//
//        int wh = Math.max(maxBinCountX, maxBinCountY);
//
//        BufferedImage image = (BufferedImage) createImage(wh, wh);
//        Graphics2D g = image.createGraphics();
//        Rectangle clipBounds = new Rectangle(0, 0, wh, wh);
//        renderer.render(0, 0, maxBinCountX, maxBinCountY, zd, displayOption, g);
//
//        return image.getScaledInstance(tw, th, Image.SCALE_SMOOTH);

        int maxBinCountX = zd0.getxGridAxis().getBinCount();
        int maxBinCountY = zd0.getyGridAxis().getBinCount();

        int wh = Math.max(maxBinCountX, maxBinCountY);

        BufferedImage image = (BufferedImage) createImage(wh, wh);
        Graphics2D g = image.createGraphics();
        Rectangle clipBounds = new Rectangle(0, 0, wh, wh);
        renderer.render(0, 0, maxBinCountX, maxBinCountY, zd0, displayOption, g);

        return image.getScaledInstance(tw, th, Image.SCALE_SMOOTH);

    }

    /**
     * Return the specified image tile, scaled by scaleFactor
     *
     * @param tileColumn column index of tile
     * @param tileRow    row index of tile
     * @return
     */
    private ImageTile getImageTile(int tileRow, int tileColumn, MainWindow.DisplayOption displayOption) {

        String key = "_" + tileRow + "_" + tileColumn + "_ " + displayOption;
        ImageTile tile = tileCache.get(key);

        if (tile == null) {

            // Image size can be smaller than tile width when zoomed out, or near the edges.
            int maxBinCountX = hic.zd.getxGridAxis().getBinCount();
            int maxBinCountY = hic.zd.getyGridAxis().getBinCount();

            if (maxBinCountX < 0 || maxBinCountY < 0) return null;

            int imageWidth = Math.min(maxBinCountX, imageTileWidth);
            int imageHeight = Math.min(maxBinCountY, imageTileWidth);

            BufferedImage image = (BufferedImage) createImage(imageWidth, imageHeight);
            Graphics2D g2D = (Graphics2D) image.getGraphics();

            final int bx0 = tileColumn * imageTileWidth;
            final int by0 = tileRow * imageTileWidth;
            renderer.render(bx0, by0, imageWidth, imageHeight, hic.zd, displayOption, g2D);

            //           if (scaleFactor > 0.999 && scaleFactor < 1.001) {
            tile = new ImageTile(image, bx0, by0);
//            } else {
//                int scaledWidth = (int) (scaleFactor * imageWidth);
//                int scaledHeight = (int) (scaleFactor * imageHeight);
//                Image scaledImage = image.getScaledInstance(scaledWidth, scaledHeight, Image.SCALE_SMOOTH);
//
//                tile = new ImageTile(scaledImage, bx0, by0);
//            }
//            tileCache.put(key, tile);
        }
        return tile;
    }


    public void clearTileCache() {
        tileCache.clear();
    }


    static class ImageTile {
        int bLeft;
        int bTop;
        Image image;

        ImageTile(Image image, int bLeft, int py0) {
            this.bLeft = bLeft;
            this.bTop = py0;
            this.image = image;
        }
    }


    class HeatmapMouseHandler extends MouseAdapter {


        DragMode dragMode = DragMode.NONE;
        private Point lastMousePoint;

        @Override
        public void mouseEntered(MouseEvent e) {
            if (straightEdgeEnabled) {
                setCursor(Cursor.getPredefinedCursor(Cursor.CROSSHAIR_CURSOR));
            } else {
                setCursor(Cursor.getDefaultCursor());
            }
        }

        @Override
        public void mouseExited(MouseEvent e) {
            hic.setCursorX(0);
            lastCursorX = 0;
            if (straightEdgeEnabled) {
                repaint();
                mainWindow.trackPanel.repaint();
            }
        }

        @Override
        public void mousePressed(final MouseEvent e) {

            if (hic.isWholeGenome()) {
                return;
            }

            if (e.isPopupTrigger()) {
                getPopupMenu().show(HeatmapPanel.this, e.getX(), e.getY());
            } else if (e.isAltDown()) {
                dragMode = DragMode.ZOOM;
            } else {
                dragMode = DragMode.PAN;
                setCursor(mainWindow.fistCursor);
            }


            lastMousePoint = e.getPoint();

        }


        @Override
        public void mouseReleased(final MouseEvent e) {

            if (e.isPopupTrigger()) {
                getPopupMenu().show(HeatmapPanel.this, e.getX(), e.getY());

            } else if (dragMode == DragMode.ZOOM && zoomRectangle != null) {

                int binX = hic.xContext.getBinOrigin() + (int) (zoomRectangle.x / hic.xContext.getScaleFactor());
                int binY = hic.yContext.getBinOrigin() + (int) (zoomRectangle.y / hic.yContext.getScaleFactor());
                int wBins = (int) (zoomRectangle.width / hic.xContext.getScaleFactor());
                int hBins = (int) (zoomRectangle.height / hic.xContext.getScaleFactor());

                int xBP0 = hic.zd.getxGridAxis().getGenomicStart(binX);
                int xBP1 = hic.zd.getxGridAxis().getGenomicEnd(binX + wBins);
                double wBP = xBP1 - xBP0;

                int yBP0 = hic.zd.getyGridAxis().getGenomicEnd(binY);
                int yBP1 = hic.zd.getyGridAxis().getGenomicEnd(binY + hBins);
                double hBP = yBP1 - yBP0;

                double newXScale = wBP / getWidth();
                double newYScale = hBP / getHeight();
                double newScale = Math.max(newXScale, newYScale);

                hic.zoomTo(xBP0, yBP0, xBP1, yBP1, newScale);

            }

            dragMode = DragMode.NONE;
            lastMousePoint = null;
            zoomRectangle = null;
            setCursor(straightEdgeEnabled ? Cursor.getPredefinedCursor(Cursor.CROSSHAIR_CURSOR) : Cursor.getDefaultCursor());
            //repaint();
        }


        @Override
        final public void mouseDragged(final MouseEvent e) {

            if (hic.zd == null || hic.isWholeGenome()) {
                return;
            }

            if (lastMousePoint == null) {
                lastMousePoint = e.getPoint();
                return;
            }

            int deltaX = e.getX() - lastMousePoint.x;
            int deltaY = e.getY() - lastMousePoint.y;
            switch (dragMode) {
                case ZOOM:
                    Rectangle lastRectangle = zoomRectangle;

                    if (deltaX == 0 || deltaY == 0) {
                        return;
                    }

                    // Constrain aspect ratio of zoom rectangle to that of panel
                    double aspectRatio = (double) getWidth() / getHeight();
                    if (deltaX * aspectRatio > deltaY) {
                        deltaY = (int) (deltaX / aspectRatio);
                    } else {
                        deltaX = (int) (deltaY * aspectRatio);
                    }


                    int x = deltaX > 0 ? lastMousePoint.x : lastMousePoint.x + (int) deltaX;
                    int y = deltaY > 0 ? lastMousePoint.y : lastMousePoint.y + (int) deltaY;
                    zoomRectangle = new Rectangle(x, y, (int) Math.abs(deltaX), (int) Math.abs(deltaY));

                    Rectangle damageRect = lastRectangle == null ? zoomRectangle : zoomRectangle.union(lastRectangle);
                    damageRect.x--;
                    damageRect.y--;
                    damageRect.width += 2;
                    damageRect.height += 2;
                    paintImmediately(damageRect);

                    break;
                default:

                    // int dx = (int) (deltaX * hic.xContext.getScale());
                    // int dy = (int) (deltaY * hic.yContext.getScale());
                    lastMousePoint = e.getPoint();    // Always save the last Point

                    hic.moveBy(-deltaX, -deltaY);

            }

        }

        @Override
        public void mouseClicked(MouseEvent e) {

            if (hic == null) return;

            if (!e.isPopupTrigger() && e.getButton() == MouseEvent.BUTTON1 && !e.isControlDown()) {

                if (hic.isWholeGenome()) {
                    int binX = hic.xContext.getBinOrigin() + (int) (e.getX() / hic.xContext.getScaleFactor());
                    int binY = hic.yContext.getBinOrigin() + (int) (e.getY() / hic.yContext.getScaleFactor());

                    int xGenome = hic.zd.getxGridAxis().getGenomicMid(binX);
                    int yGenome = hic.zd.getyGridAxis().getGenomicMid(binY);

                    Chromosome xChrom = null;
                    Chromosome yChrom = null;
                    for (int i = 0; i < chromosomeBoundaries.length; i++) {
                        if (xChrom == null && chromosomeBoundaries[i] > xGenome) {
                            xChrom = hic.getChromosomes()[i + 1];
                        }
                        if (yChrom == null && chromosomeBoundaries[i] > yGenome) {
                            yChrom = hic.getChromosomes()[i + 1];
                        }
                        if (xChrom != null && yChrom != null) {
                            mainWindow.setSelectedChromosomes(xChrom, yChrom);
                        }
                    }


                } else if (e.getClickCount() > 1) {
                    // Double click,  zoom and center on click location
                    int currentZoom = hic.xContext.getZoom();
                    final int newZoom = e.isAltDown()
                            ? Math.max(currentZoom - 1, 1)
                            : Math.min(11, currentZoom + 1);


                    int centerBinX = hic.xContext.getBinOrigin() + (int) (e.getX() / hic.xContext.getScaleFactor());
                    int centerBinY = hic.yContext.getBinOrigin() + (int) (e.getY() / hic.yContext.getScaleFactor());

                    int xGenome = hic.zd.getxGridAxis().getGenomicMid(centerBinX);
                    int yGenome = hic.zd.getyGridAxis().getGenomicMid(centerBinY);

                    hic.setZoom(newZoom, xGenome, yGenome);

                } else {

                    if (hic == null) return;

                    int centerBinX = hic.xContext.getBinOrigin() + (int) (e.getX() / hic.xContext.getScaleFactor());
                    int centerBinY = hic.yContext.getBinOrigin() + (int) (e.getY() / hic.yContext.getScaleFactor());

                    int xGenome = hic.zd.getxGridAxis().getGenomicMid(centerBinX);
                    int yGenome = hic.zd.getyGridAxis().getGenomicMid(centerBinY);

                    System.out.println("Bin position:  " + centerBinX + "\t" + centerBinY);
                    System.out.println("Genome position:  " + xGenome + "\t" + yGenome);


                    //If IGV is running open on loci
                    if (e.isShiftDown()) {

//                        String chr1 = hic.xContext.getChromosome().getName();
//                        int leftX = (int) hic.xContext.getChromosomePosition(0);
//                        int wX = (int) (hic.xContext.getScale() * getWidth());
//                        int rightX = leftX + wX;
//
//                        String chr2 = hic.yContext.getChromosome().getName();
//                        int leftY = (int) hic.yContext.getChromosomePosition(0);
//                        int wY = (int) (hic.xContext.getScale() * getHeight());
//                        int rightY = leftY + wY;
//
//                        String locus1 = "chr" + chr1 + ":" + leftX + "-" + rightX;
//                        String locus2 = "chr" + chr2 + ":" + leftY + "-" + rightY;
//
//                        IGVUtils.sendToIGV(locus1, locus2);
                    }
                }

            }
        }

        @Override
        public void mouseMoved(MouseEvent e) {
            if (hic.xContext != null && hic.zd != null) {


                if (straightEdgeEnabled) {
                    hic.setCursorX(e.getX());
                    Rectangle damageRect = new Rectangle();
                    damageRect.x = (lastCursorX > 0 ? Math.min(lastCursorX, e.getX()) : e.getX()) - 1;
                    damageRect.y = 0;
                    damageRect.width = lastCursorX > 0 ? Math.abs(e.getX() - lastCursorX) + 2 : 2;
                    damageRect.height = getHeight();

                    paintImmediately(damageRect);

                    damageRect.height = mainWindow.trackPanel.getHeight();
                    mainWindow.trackPanel.paintImmediately(damageRect);

                    lastCursorX = e.getX();

                }


                HiCGridAxis xGridAxis = hic.zd.getxGridAxis();
                HiCGridAxis yGridAxis = hic.zd.getyGridAxis();
                if (hic.isWholeGenome()) {

                    int binX = hic.xContext.getBinOrigin() + (int) (e.getX() / hic.xContext.getScaleFactor());
                    int binY = hic.xContext.getBinOrigin() + (int) (e.getY() / hic.xContext.getScaleFactor());

                    int xGenome = xGridAxis.getGenomicStart(binX);
                    int yGenome = yGridAxis.getGenomicStart(binY);

                    Chromosome xChrom = null;
                    Chromosome yChrom = null;
                    for (int i = 0; i < chromosomeBoundaries.length; i++) {
                        if (xChrom == null && chromosomeBoundaries[i] > xGenome) {
                            xChrom = hic.getChromosomes()[i + 1];
                        }
                        if (yChrom == null && chromosomeBoundaries[i] > yGenome) {
                            yChrom = hic.getChromosomes()[i + 1];
                        }
                        if (xChrom != null && yChrom != null) {

                            int leftBoundaryX = xChrom.getIndex() == 1 ? 0 : chromosomeBoundaries[xChrom.getIndex() - 2];
                            int leftBoundaryY = yChrom.getIndex() == 1 ? 0 : chromosomeBoundaries[yChrom.getIndex() - 2];


                            int xChromPos = (int) ((xGenome - leftBoundaryX) * 1000);
                            int yChromPos = (int) ((yGenome - leftBoundaryY) * 1000);

                            StringBuffer txt = new StringBuffer();
                            txt.append(xChrom.getName());
                            txt.append(":");
                            txt.append(String.valueOf(xChromPos));
                            txt.append("<html>");
                            txt.append("<br>");

                            txt.append(yChrom.getName());
                            txt.append(":");
                            txt.append(String.valueOf(yChromPos));

                            setToolTipText(txt.toString());
                            return;

                        }
                    }

                } else {
                    int binX = hic.xContext.getBinOrigin() + (int) (e.getX() / hic.xContext.getScaleFactor());
                    int binY = hic.xContext.getBinOrigin() + (int) (e.getY() / hic.xContext.getScaleFactor());

                    int xGenome = xGridAxis.getGenomicStart(binX);
                    int yGenome = yGridAxis.getGenomicStart(binY);
                    int xGenomeEnd = xGridAxis.getGenomicEnd(binX);
                    int yGenomeEnd = yGridAxis.getGenomicEnd(binY);

                    //int binX = (int) ((mainWindow.xContext.getOrigin() + e.getX() * mainWindow.xContext.getScale()) / getBinWidth());
                    //int binY = (int) ((mainWindow.yContext.getOrigin() + e.getY() * mainWindow.yContext.getScale()) / getBinWidth());
                    StringBuffer txt = new StringBuffer();
                    txt.append("<html>");
                    txt.append(hic.xContext.getChromosome().getName());
                    txt.append(":");
                    txt.append(formatter.format(xGenome));
                    txt.append("-");
                    txt.append(formatter.format(xGenomeEnd));
                    if (xGridAxis instanceof HiCFragmentAxis) {
                        txt.append("  (" + formatter.format(binX) + ")");
                    }

                    txt.append("<br>");
                    txt.append(hic.yContext.getChromosome().getName());
                    txt.append(":");
                    txt.append(formatter.format(yGenome));
                    txt.append("-");
                    txt.append(formatter.format(yGenomeEnd));

                    if (yGridAxis instanceof HiCFragmentAxis) {
                        txt.append("  (" + formatter.format(binY) + ")");
                    }

                    setToolTipText(txt.toString());
                }
            }
        }
    }


    boolean straightEdgeEnabled = false;

    JidePopupMenu getPopupMenu() {

        JidePopupMenu menu = new JidePopupMenu();

        final JCheckBoxMenuItem mi = new JCheckBoxMenuItem("Enable straight edge");
        mi.setSelected(straightEdgeEnabled);
        mi.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if (mi.isSelected()) {
                    straightEdgeEnabled = true;
                    setCursor(Cursor.getPredefinedCursor(Cursor.CROSSHAIR_CURSOR));
                } else {
                    straightEdgeEnabled = false;
                    lastCursorX = 0;
                    hic.setCursorX(0);
                    setCursor(Cursor.getDefaultCursor());
                    repaint();
                    mainWindow.trackPanel.repaint();
                }

            }
        });
        menu.add(mi);

        final JMenuItem mi2 = new JMenuItem("Goto fragment...");
        mi2.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                String fragmentString = JOptionPane.showInputDialog(HeatmapPanel.this, "Enter fragment range in the form x:y");
                if (fragmentString != null) {
                    String[] tokens = Globals.colonPattern.split(fragmentString);
                    int x = Integer.parseInt(tokens[0].replace(",", ""));
                    int y = (tokens.length > 1) ?  Integer.parseInt(tokens[1].replace(",", "")) : x;
                    hic.centerFragment(x, y);

                }
            }
        });
        if (hic != null) {
            menu.add(mi2);
        }


        return menu;

    }

}
