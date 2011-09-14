package org.broad.igv.hic;

import org.apache.batik.dom.util.HashTable;
import org.broad.igv.hic.data.MatrixZoomData;
import org.broad.igv.util.ObjectCache;
import org.jfree.ui.RectangleEdge;

import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.Serializable;
import java.util.Hashtable;

/**
 * @author jrobinso
 * @date Aug 2, 2010
 */
public class HeatmapPanel extends JComponent implements Serializable {

    private MainWindow mainWindow;
    private double maxCount = 200.0 / (1 * 1);
    private int binWidth = 2;
    private int imageTileWidth = 500;
    ObjectCache<String, ImageTile> tileCache = new ObjectCache(100);

    HeatmapRenderer renderer = new HeatmapRenderer();


    public HeatmapPanel() {

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
    }

    public Image getThumbnailImage(MatrixZoomData zd, int tw, int th) {

        int w =  Math.min(getWidth(), mainWindow.xContext.getScreenPosition(mainWindow.chr1.getSize()));
        int h = Math.min(getHeight(), mainWindow.yContext.getScreenPosition(mainWindow.chr2.getSize()));

        BufferedImage image = (BufferedImage) createImage(w, h);
        Rectangle bounds = new Rectangle(0, 0, w, h);
        Graphics g = image.createGraphics();
        renderer.render(mainWindow.xContext.getOrigin(), mainWindow.yContext.getOrigin(), mainWindow.zd,
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
}
