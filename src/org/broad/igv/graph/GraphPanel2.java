/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.graph;

import javax.swing.*;
import java.awt.*;

/**
 * @author jrobinso
 * @date Oct 12, 2010
 */
public class GraphPanel2 extends JPanel {

    private Graph graph;
    int topMargin = 50;
    int leftMargin = 20;

    int nodeWidth = 5;
    int nodeHeight = 5;
    int ySpacing = 20;


    public void setGraph(Graph graph) {
        this.graph = graph;
    }

    @Override
    protected void paintComponent(Graphics graphics) {

        Graphics2D g = (Graphics2D) graphics;

        ((Graphics2D) g).setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        // Just paint nodes (no edges)
        if (graph != null) {

            // Set X scale to fit graph in panel
            float dx = graph.getMaxX() - graph.getMinX();
            float scale = (getWidth() - 2 * leftMargin) / dx;

            for (Node node : graph.getNodes()) {

                int pixelX = getXPosition(scale, node.getCellX());
                int pixelY = getYPosition(node);

                g.setColor(Color.BLUE);
                g.fillRect(pixelX, pixelY, nodeWidth, nodeHeight);

            }

            // Paint edges
            for (SubGraph sg : graph.getSubGraphs()) {
                for (Edge edge : sg.getEdges()) {

                    // Parent right position
                    Node parent = edge.getParent();
                    int px1 = getXPosition(scale, parent.getCellX()) + nodeWidth;
                    int py1 = getYPosition(parent) + nodeHeight / 2;

                    // Child left position
                    Node child = edge.getChild();
                    int px2 = getXPosition(scale, child.getCellX());
                    int py2 = getYPosition(child) + nodeHeight / 2;

                    g.setColor(edge.getColor());
                    g.drawLine(px1, py1, px2, py2);

                    g.setColor(Color.red);
                    g.fillRect(px2 - 1, py2 - 1, 2, 2);


                }
            }
        }


    }

    private int getYPosition(Node node) {
        return topMargin + node.getCellY() * ySpacing;
    }

    private int getXPosition(float scale, int x) {
        return leftMargin + (int) (scale * (x - graph.getMinX()));
    }


}
