package org.broad.igv.hic;

import org.broad.igv.hic.data.MatrixZoomData;
import org.broad.igv.util.ObjectCache;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.io.Serializable;

/**
 * @author jrobinso
 * @date Aug 2, 2010
 */
public class HeatmapPanel extends JComponent implements Serializable {

    private MainWindow mainWindow;
    private double maxCount = 20000.0 / (1 * 1);
    private int imageTileWidth = 500;
    ObjectCache<String, ImageTile> tileCache = new ObjectCache(100);
    private Rectangle zoomRectangle;

    HeatmapRenderer renderer = new HeatmapRenderer();


    public HeatmapPanel() {
        final HeatmapMouseHandler mouseHandler = new HeatmapMouseHandler();
        addMouseListener(mouseHandler);
        addMouseMotionListener(mouseHandler);
    }

    @Override
    protected void paintComponent(Graphics g) {

        if (mainWindow != null && mainWindow.zd != null) {

            // 1  get image in "bin" coordinates.

            // bp / bin
            int binSize = mainWindow.zd.getBinSize();

            int originX = mainWindow.xContext.getOrigin();
            int originY = mainWindow.yContext.getOrigin();

            int bLeft = (int) (originX / binSize);
            int bRight = bLeft + (int) ((getWidth() * mainWindow.xContext.getScale()) / binSize);
            int bTop = (int) (originY / binSize);
            int bBottom = bTop + (int) ((getWidth() * mainWindow.xContext.getScale()) / binSize);

            // tile coordinates
            int tLeft = bLeft / imageTileWidth;
            int tRight = bRight / imageTileWidth;
            int tTop = bTop / imageTileWidth;
            int tBottom = bBottom / imageTileWidth;

            // scale factor -- for now we are using the same scale for x & y (square pixels)
            double pixelsPerBin = binSize / mainWindow.xContext.getScale();

            for (int i = tLeft; i <= tRight; i++) {
                for (int j = tTop; j <= tBottom; j++) {

                    ImageTile tile = getImageTile(i, j, pixelsPerBin);

                    int pxOffset = (int) ((tile.bLeft - bLeft) * pixelsPerBin);
                    int pyOffset = (int) ((tile.bTop - bTop) * pixelsPerBin);

                    g.drawImage(tile.image, pxOffset, pyOffset, null);

                }
            }

        }

        if (zoomRectangle != null) {
            ((Graphics2D) g).draw(zoomRectangle);
        }
    }

    public Image getThumbnailImage(MatrixZoomData zd, int tw, int th) {

        int w = Math.min(getWidth(), mainWindow.xContext.getScreenPosition(mainWindow.xContext.getChromosome().getSize()));
        int h = Math.min(getHeight(), mainWindow.yContext.getScreenPosition(mainWindow.yContext.getChromosome().getSize()));
        int wh = Math.max(w, h);

        BufferedImage image = (BufferedImage) createImage(wh, wh);
        Rectangle bounds = new Rectangle(0, 0, w, h);
        Graphics g = image.createGraphics();

        // TODO -- get # of bins
        int nBins = 500;

        renderer.render(0, 0, 500, 500, zd, maxCount, g, getBackground());

        return image.getScaledInstance(tw, th, Image.SCALE_SMOOTH);

    }

    private ImageTile getImageTile(int i, int j, double scaleFactor) {
        String key = "_" + i + "_" + j;
        ImageTile tile = tileCache.get(key);
        if (tile == null) {

            // Image size can be smaller than tile width when zoomed out
            int maxBinCountX = mainWindow.xContext.getChrLength() / mainWindow.zd.getBinSize() + 1;
            int maxBinCountY = mainWindow.yContext.getChrLength() / mainWindow.zd.getBinSize() + 1;

            int imageWidth = Math.min(maxBinCountX, imageTileWidth);
            int imageHeight = Math.min(maxBinCountY, imageTileWidth);


            BufferedImage image = (BufferedImage) createImage(imageWidth, imageHeight);
            Graphics2D g2D = (Graphics2D) image.getGraphics();

            final int bx0 = i * imageTileWidth;
            final int by0 = j * imageTileWidth;
            renderer.render(bx0, by0, imageWidth, imageHeight, mainWindow.zd, maxCount, g2D, getBackground());

            if (scaleFactor > 0.999 && scaleFactor < 1.001) {
                tile = new ImageTile(image, bx0, by0);
            } else {
                int scaledWidth = (int) (scaleFactor * imageWidth);
                int scaledHeight = (int) (scaleFactor * imageHeight);
                Image scaledImage = image.getScaledInstance(scaledWidth, scaledHeight, Image.SCALE_SMOOTH);

                tile = new ImageTile(scaledImage, bx0, by0);
            }
            tileCache.put(key, tile);
        }
        return tile;
    }

    public void setMaxCount(double maxCount) {
        this.maxCount = maxCount;
        tileCache.clear();
    }

    public MainWindow getMainWindow() {
        return mainWindow;
    }

    public void setMainWindow(MainWindow mainWindow) {
        this.mainWindow = mainWindow;
    }


    @Override
    public void setSize(int i, int i1) {
        super.setSize(i, i1);    //To change body of overridden methods use File | Settings | File Templates.
    }

    @Override
    public void setSize(Dimension dimension) {
        super.setSize(dimension);    //To change body of overridden methods use File | Settings | File Templates.
    }

    @Override
    public void setBounds(Rectangle rectangle) {
        super.setBounds(rectangle);    //To change body of overridden methods use File | Settings | File Templates.
    }

    @Override
    public void setBounds(int i, int i1, int i2, int i3) {
        super.setBounds(i, i1, i2, i3);    //To change body of overridden methods use File | Settings | File Templates.
    }

    public void clearTileCache() {
        System.out.println("Clear cache");
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

        int dragMode = 0;  // 0 => none, 1 => pan,  2 => zoom
        private Point lastMousePoint;


        @Override
        public void mousePressed(final MouseEvent e) {

            if (e.isAltDown()) {
                dragMode = 2;
            } else {
                dragMode = 1;
                setCursor(mainWindow.fistCursor);
            }

            lastMousePoint = e.getPoint();

        }


        @Override
        public void mouseReleased(final MouseEvent e) {

            if (dragMode == 2) {

                double xBP = mainWindow.xContext.getOrigin() + zoomRectangle.getX() * mainWindow.xContext.getScale();
                double yBP = mainWindow.yContext.getOrigin() + zoomRectangle.getY() * mainWindow.yContext.getScale();
                double wBP = zoomRectangle.width * mainWindow.xContext.getScale();
                double hBP = zoomRectangle.height * mainWindow.yContext.getScale();

                double newXScale = wBP / getWidth();
                double newYScale = hBP / getHeight();
                double newScale = Math.max(newXScale, newYScale);

                mainWindow.zoomTo(xBP, yBP, newScale);

            }

            dragMode = 0;
            lastMousePoint = null;
            zoomRectangle = null;
            setCursor(Cursor.getDefaultCursor());
            repaint();
        }


        @Override
        final public void mouseDragged(final MouseEvent e) {

            if (mainWindow.zd == null) {
                return;
            }

            if (lastMousePoint == null) {
                lastMousePoint = e.getPoint();
                return;
            }

            double deltaX = e.getX() - lastMousePoint.x;
            double deltaY = e.getY() - lastMousePoint.y;
            switch (dragMode) {
                case 2:

                    Rectangle lastRectangle = zoomRectangle;

                    if (deltaX == 0 || deltaY == 0) {
                        return;
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

                    // Size i
                    int dx = (int) (deltaX * mainWindow.xContext.getScale());
                    int dy = (int) (deltaY * mainWindow.yContext.getScale());
                    lastMousePoint = e.getPoint();    // Always save the last Point

                    mainWindow.moveBy(-dx, -dy);

            }

        }

        @Override
        public void mouseClicked(MouseEvent e) {

            if (!e.isPopupTrigger() && (e.getClickCount() > 1)) {
                int currentZoom = mainWindow.xContext.getZoom();
                final int newZoom = e.isAltDown()
                        ? Math.max(currentZoom - 1, 1)
                        : Math.min(11, currentZoom + 1);

                int centerLocationX = mainWindow.xContext.getOrigin() + (int) (e.getX() * mainWindow.xContext.getScale());
                int centerLocationY = mainWindow.yContext.getOrigin() + (int) (e.getY() * mainWindow.yContext.getScale());

                mainWindow.setZoom(newZoom, centerLocationX, centerLocationY);

            }
        }

        @Override
        public void mouseMoved(MouseEvent e) {
            if (mainWindow.xContext != null && mainWindow.zd != null) {
                Rectangle visRect = getVisibleRect();
                final int binSize = mainWindow.zd.getBinSize();

                int posX = (int) (e.getX() / mainWindow.xContext.getScale());
                int posY = (int) (e.getY() / mainWindow.yContext.getScale());

                //int binX = (int) ((mainWindow.xContext.getOrigin() + e.getX() * mainWindow.xContext.getScale()) / getBinWidth());
                //int binY = (int) ((mainWindow.yContext.getOrigin() + e.getY() * mainWindow.yContext.getScale()) / getBinWidth());
                StringBuffer txt = new StringBuffer();
                txt.append("<html>");
                txt.append(mainWindow.xContext.getChromosome().getName());
                txt.append(":");
                txt.append(String.valueOf(posX));
                txt.append("<br>");
                txt.append(mainWindow.yContext.getChromosome().getName());
                txt.append(":");
                txt.append(String.valueOf(posY));
                setToolTipText(txt.toString());
            }
        }
    }

}
