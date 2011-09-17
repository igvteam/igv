package org.broad.igv.hic;

import org.apache.batik.dom.util.HashTable;
import org.broad.igv.hic.data.MatrixZoomData;
import org.broad.igv.util.ObjectCache;
import org.jfree.ui.RectangleEdge;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.io.Serializable;
import java.util.Hashtable;

/**
 * @author jrobinso
 * @date Aug 2, 2010
 */
public class HeatmapPanel extends JComponent implements Serializable {

    private MainWindow mainWindow;
    private double maxCount = 20000.0 / (1 * 1);
    private int binWidth = 2;
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

            int originX = mainWindow.xContext.getOrigin();
            int originY = mainWindow.yContext.getOrigin();

            int px0 = (int) (originX / mainWindow.xContext.getScale());
            int px1 = px0 + getWidth();

            int py0 = (int) (originY / mainWindow.yContext.getScale());
            int py1 = py0 + getHeight();

            int tx0 = px0 / imageTileWidth;
            int tx1 = px1 / imageTileWidth;

            int ty0 = py0 / imageTileWidth;
            int ty1 = py1 / imageTileWidth;

            int zoomLevel = mainWindow.zd.getZoom();
            for (int i = tx0; i <= tx1; i++) {
                for (int j = ty0; j <= ty1; j++) {
                    ImageTile tile = getImageTile(zoomLevel, i, j);

                    int pxOffset = tile.px0 - px0;
                    int pyOffset = tile.py0 - py0;

                    g.drawImage(tile.image, pxOffset, pyOffset, null);

                }
            }

        }

        if(zoomRectangle != null) {
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
        renderer.render(mainWindow.xContext.getOrigin(), mainWindow.yContext.getOrigin(), zd,
                binWidth, maxCount, g, bounds, getBackground());

        return image.getScaledInstance(tw, th, Image.SCALE_SMOOTH);

    }


    private ImageTile getImageTile(int z, int i, int j) {
        String key = z + "_" + i + "_" + j;
        ImageTile tile = tileCache.get(key);
        if (tile == null) {

            BufferedImage image = (BufferedImage) createImage(imageTileWidth, imageTileWidth);
            Graphics2D g2D = (Graphics2D) image.getGraphics();

            final int px0 = i * imageTileWidth;
            int x0 = (int) (px0 * mainWindow.xContext.getScale());
            final int py0 = j * imageTileWidth;
            int y0 = (int) (py0 * mainWindow.yContext.getScale());
            Rectangle bounds = new Rectangle(0, 0, imageTileWidth, imageTileWidth);
            renderer.render(x0, y0, mainWindow.zd, binWidth, maxCount, g2D, bounds, getBackground());

            tile = new ImageTile(image, px0, py0);
            tileCache.put(key, tile);
        }
        return tile;
    }

    public void setMaxCount(double maxCount) {
        this.maxCount = maxCount;
        tileCache.clear();
    }

    public int getBinWidth() {
        return binWidth;
    }

    public MainWindow getMainWindow() {
        return mainWindow;
    }

    public void setMainWindow(MainWindow mainWindow) {
        this.mainWindow = mainWindow;
    }

    public void setBinWidth(int binWidth) {
        System.out.println("bin width = " + binWidth);
        this.binWidth = binWidth;
        tileCache.clear();
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
        tileCache.clear();
    }


    static class ImageTile {
        int px0;
        int py0;
        BufferedImage image;

        ImageTile(BufferedImage image, int px0, int py0) {
            this.px0 = px0;
            this.py0 = py0;
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
            dragMode = 0;
            lastMousePoint = null;
            zoomRectangle = null;
            setCursor(Cursor.getDefaultCursor());
            repaint();
        }


        @Override
        final public void mouseDragged(final MouseEvent e) {

            if( mainWindow.zd == null) {
                return;
            }

            if (lastMousePoint == null) {
                lastMousePoint = e.getPoint();
                return;
            }

            double deltaX =  e.getX() - lastMousePoint.x;
            double deltaY =  e.getY() - lastMousePoint.y;
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
                    int dx = (int) (deltaX * mainWindow.zd.getBinSize() / getBinWidth());
                    int dy = (int) (deltaY * mainWindow.zd.getBinSize() / getBinWidth());
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

                int centerLocationX = mainWindow.xContext.getOrigin() + e.getX() * mainWindow.zd.getBinSize() / getBinWidth();
                int centerLocationY = mainWindow.yContext.getOrigin() + e.getY() * mainWindow.zd.getBinSize() / getBinWidth();

                mainWindow.setZoom(newZoom, centerLocationX, centerLocationY);

            }
        }

        @Override
        public void mouseMoved(MouseEvent e) {
            if (mainWindow.xContext != null && mainWindow.zd != null) {
                Rectangle visRect = getVisibleRect();
                final int binSize = mainWindow.zd.getBinSize();
                int binX = (mainWindow.xContext.getOrigin() / binSize) + (e.getX() - visRect.x) / getBinWidth();
                int binY = (mainWindow.yContext.getOrigin() / binSize) + (e.getY() - visRect.y) / getBinWidth();
                StringBuffer txt = new StringBuffer();
                txt.append("<html>");
                txt.append(mainWindow.xContext.getChromosome().getName());
                txt.append(":");
                txt.append(String.valueOf((binX - 1) * binSize));
                txt.append("-");
                txt.append(String.valueOf(binX * binSize));
                txt.append("<br>");
                txt.append(mainWindow.yContext.getChromosome().getName());
                txt.append(":");
                txt.append(String.valueOf((binY - 1) * binSize));
                txt.append("-");
                txt.append(String.valueOf(binY * binSize));
                setToolTipText(txt.toString());
            }
        }
    }

}
