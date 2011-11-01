package org.broad.igv.charts;

import org.broad.igv.util.stats.KMPlotFrame;
import org.omg.CORBA.Bounds;

import java.awt.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Simple geometric class to support lasso selection
 *
 * @author Jim Robinson
 * @date 10/31/11
 */
public class SelectionPath {

    List<Point>  points;
    Rectangle bounds;
    int maxX;
    int maxY;

    public SelectionPath() {
        points = new ArrayList<Point>(100);
    }

    public void addPoint(Point p) {

        if(bounds == null) {
            // Be generous with bounds, its only use is to set clip rects for repaints
            bounds = new Rectangle(p.x-10, p.y-10, 20, 20);
            maxX = p.x + 10;
            maxY = p.y + 10;
        }
        else {
            maxX = Math.max(maxX, p.x + 10);
            bounds.x = Math.min(bounds.x, p.x - 10);
            bounds.width = maxX - bounds.x;

            maxY = Math.max(maxY, p.y + 10);
            bounds.y = Math.min(bounds.y, p.y - 10);
            bounds.height = maxY - bounds.y;
        }
        points.add(p);
    }

    public List<Point> getPoints() {
        return points;
    }

    public Rectangle getBounds() {
        return bounds;
    }

    public int size() {
        return points.size();
    }
}
