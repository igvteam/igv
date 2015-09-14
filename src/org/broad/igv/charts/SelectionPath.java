/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

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
