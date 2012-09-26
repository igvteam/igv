package org.broad.igv.hic;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.awt.geom.*;
import java.io.Serializable;

/**
 * @author jrobinso
 * @date Aug 2, 2010
 */
public class ThumbnailPanel extends JComponent implements Serializable {


    private MainWindow mainWindow;
    private HiC hic;

    Image image;

    public static final AlphaComposite ALPHA_COMP = AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.75f);
    private double yScale;
    private double xScale;
    Point lastPoint = null;

    private Rectangle innerRectangle;


    public ThumbnailPanel(MainWindow mainWindow, HiC model) {

        this.mainWindow = mainWindow;
        this.hic = model;

        addMouseListener(new MouseAdapter() {


            @Override
            public void mouseClicked(MouseEvent mouseEvent) {
                int xBP = (int) (mouseEvent.getX() * xScale);
                int yBP = (int) (mouseEvent.getY() * yScale);
                hic.center(xBP, yBP);
            }

            @Override
            public void mousePressed(MouseEvent mouseEvent) {
                if (innerRectangle != null && innerRectangle.contains(mouseEvent.getPoint())) {
                    lastPoint = mouseEvent.getPoint();
                    setCursor(MainWindow.fistCursor);
                }

            }

            @Override
            public void mouseReleased(MouseEvent mouseEvent) {
                lastPoint = null;
                setCursor(Cursor.getDefaultCursor());
            }
        });

        addMouseMotionListener(new MouseMotionAdapter() {
            @Override
            public void mouseDragged(MouseEvent mouseEvent) {
                if (lastPoint != null) {
                    int dxBP = ((int) ((mouseEvent.getX() - lastPoint.x) * xScale));
                    int dyBP = ((int) ((mouseEvent.getY() - lastPoint.y) * xScale));
                    hic.moveBy(dxBP, dyBP);
                    lastPoint = mouseEvent.getPoint();
                }


            }
        });

    }

    public void setImage(Image image) {
        this.image = image;
        xScale = ((double) hic.zd.getxGridAxis().getBinCount()) / getWidth();
        yScale = ((double) hic.zd.getyGridAxis().getBinCount()) / getWidth();
        ;
    }

    public String getName() {
        return null;
    }

    public void setName(String nm) {

    }

    @Override
    protected void paintComponent(Graphics g) {

        if (image != null) {
            g.drawImage(image, 0, 0, null);
            renderVisibleWindow((Graphics2D) g);
        }
    }

    private void renderVisibleWindow(Graphics2D g) {


        if (hic != null && hic.xContext != null) {

            Rectangle outerRectangle = new Rectangle(0, 0, getBounds().width, getBounds().height);

            int wPixels = mainWindow.getHeatmapPanel().getWidth();
            int hPixels = mainWindow.getHeatmapPanel().getHeight();

            int originX = hic.xContext.getBinOrigin();
            int x = (int) (originX / xScale);

            int originY = hic.yContext.getBinOrigin();
            int y = (int) (originY / yScale);

            int wBins = wPixels;  // TODO -- scale @ lower resolution
            int w = (int) (wBins / xScale);

            int yBins = hPixels;
            int h = (int) (yBins / yScale);

            if (w < 4) {
                int delta = 4 - w;
                x -= delta / 2;
                w = 4;
            }
            if (h < 4) {
                int delta = 4 - h;
                y -= delta / 2;
                h = 4;
            }

            innerRectangle = new Rectangle(x, y, w, h);
            Shape shape = new SquareDonut(outerRectangle, innerRectangle);

            g.setColor(Color.gray);
            g.setComposite(ALPHA_COMP);
            g.fill(shape);

            g.draw(innerRectangle);
        }
    }

    private static class SquareDonut implements Shape {
        private final Area area;

        public SquareDonut(Rectangle outerRectangle, Rectangle innerRectangle) {
            this.area = new Area(outerRectangle);
            area.subtract(new Area(innerRectangle));
        }

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
    }
}
