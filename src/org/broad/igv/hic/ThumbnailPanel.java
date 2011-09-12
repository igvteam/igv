package org.broad.igv.hic;

import org.broad.igv.hic.data.*;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.*;
import java.io.Serializable;

/**
 * @author jrobinso
 * @date Aug 2, 2010
 */
public class ThumbnailPanel extends JComponent implements Serializable {

    private MainWindow mainWindow;
    private MatrixZoomData zd;
    private int maxCount = 500;
    private int binWidth = 1;
    private Context context;

    HeatmapRenderer renderer = new HeatmapRenderer();


    public ThumbnailPanel() {

    }

    public void setMainWindow(MainWindow mainWindow) {
        this.mainWindow = mainWindow;
    }

    public String getName() {
        return null;
    }

    public void setName(String nm) {

    }


    @Override
    protected void paintComponent(Graphics g) {

        if (getZd() != null) {
            Rectangle bounds = this.getVisibleRect();

            renderer.render(0, 0, getZd(), getBinWidth(), getMaxCount(), g, bounds, getBackground());
            renderVisibleWindow((Graphics2D) g);
        }
    }

    private void renderVisibleWindow(Graphics2D g) {
        if (mainWindow != null && mainWindow.xContext != null && mainWindow.xContext.getVisibleWidth() > 0) {

            int bw = getBinWidth();
            int nBins = mainWindow.getLen() / zd.getBinSize();
            int effectiveWidth = nBins * bw;

            int w = (int) ((((double) mainWindow.xContext.getVisibleWidth()) / mainWindow.getLen()) * effectiveWidth);
            int h = (int) ((((double) mainWindow.yContext.getVisibleWidth()) / mainWindow.getLen()) * effectiveWidth);
            int x = (int) ((((double) mainWindow.xContext.getOrigin()) / mainWindow.getLen()) * effectiveWidth);
            int y = (int) ((((double) mainWindow.yContext.getOrigin()) / mainWindow.getLen()) * effectiveWidth);

            Rectangle outerRectangle = new Rectangle(0, 0, getBounds().width, getBounds().height);
            Rectangle innerRectangle = new Rectangle(x, y, w, h);
            final Area area = new Area(outerRectangle);
            area.subtract(new Area(innerRectangle));

            Shape shape = new Shape() {
                public Rectangle getBounds() {
                    return area.getBounds();
                }

                public Rectangle2D getBounds2D() {
                    return area.getBounds2D();
                }

                public boolean contains(double v, double v1) {
                    return area.contains(v, v1);
                }

                public boolean contains(Point2D point2D) {
                    return area.contains(point2D);
                }

                public boolean intersects(double v, double v1, double v2, double v3) {
                    return area.intersects(v, v1, v2, v3);
                }

                public boolean intersects(Rectangle2D rectangle2D) {
                    return area.intersects(rectangle2D);
                }

                public boolean contains(double v, double v1, double v2, double v3) {
                    return area.contains(v, v1, v2, v3);
                }

                public boolean contains(Rectangle2D rectangle2D) {
                    return area.contains(rectangle2D);
                }

                public PathIterator getPathIterator(AffineTransform affineTransform) {
                    return area.getPathIterator(affineTransform);
                }

                public PathIterator getPathIterator(AffineTransform affineTransform, double v) {
                    return area.getPathIterator(affineTransform, v);
                }
            };

            g.setColor(Color.gray);
            AlphaComposite alphaComp = AlphaComposite.getInstance(
                    AlphaComposite.SRC_OVER, 0.75f);
            g.setComposite(alphaComp);
            g.fill(shape);

            //g.setColor(Color.black);
            g.draw(innerRectangle);
        }
    }

    public MatrixZoomData getZd() {
        return zd;
    }

    public void setZd(MatrixZoomData zd) {
        this.zd = zd;
        // TODO -- draw image here
    }

    public int getMaxCount() {
        return maxCount;
    }

    public void setMaxCount(int maxCount) {
        this.maxCount = maxCount;
    }

    public int getBinWidth() {
        return binWidth;
    }

    public void setBinWidth(int binWidth) {
        this.binWidth = binWidth;
    }

    public Context getContext() {
        return context;
    }

    public void setContext(Context context) {
        this.context = context;
    }
}
