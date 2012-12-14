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


/*
* LocationPanel.java
*
* Created on September 11, 2007, 2:29 PM
*
* To change this template, choose Tools | Template Manager
* and open the template in the editor.
*/
package org.broad.igv.ui.panel;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.Globals;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.Cytoband;
import org.broad.igv.renderer.CytobandRenderer;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.WaitCursorManager;

import javax.swing.*;
import javax.swing.event.MouseInputAdapter;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 */
public class CytobandPanel extends JPanel {

    private static int fontHeight = 10;
    private static int bandHeight = 10;
    private static String fontFamilyName = "Lucida Sans";
    private boolean isDragging = false;

    private double viewOrigin;
    private double viewEnd;

    double cytobandScale;
    ReferenceFrame frame;
    private Rectangle currentRegionRect;
    private CytobandRenderer cytobandRenderer;
    private List<Cytoband> currentCytobands;

    public CytobandPanel(ReferenceFrame frame) {
        this(frame, true);
    }


    public CytobandPanel(ReferenceFrame frame, boolean mouseable) {

        this.frame = frame;
        viewOrigin = frame.getOrigin();
        viewEnd = frame.getEnd();

        FontManager.getFont(fontHeight);
        setFont(new Font(fontFamilyName, Font.BOLD, fontHeight));
        if (mouseable) {
            initMouseAdapter();
        }
        cytobandRenderer = (new CytobandRenderer());
    }

    
    @Override
    protected void paintComponent(Graphics g) {

        super.paintComponent(g);

        ((Graphics2D) g).setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        if (frame.getChrName().equals(Globals.CHR_ALL) || getWidth() < 10) {
            //Graphics g2 = g.create();
            //g2.setFont(FontManager.getScalableFont(Font.ITALIC, 12));
            //String text = "Whole genome view.  To jump to a chromosome click on its label.";
            //g2.drawString(text, 20, getHeight() - 5);
            return;
        }

        int dataPanelWidth = frame.getWidthInPixels();
        Rectangle cytoRect = new Rectangle(0, 10, dataPanelWidth, bandHeight);

        Chromosome chromosome = getReferenceFrame().getChromosome();
        if (chromosome == null) {
            return;
        }
        currentCytobands = chromosome.getCytobands();
        if (currentCytobands == null) {
            return;
        }

        cytobandRenderer.draw(currentCytobands, g, cytoRect, frame);

        int chromosomeLength = getReferenceFrame().getChromosomeLength();
        cytobandScale = ((double) chromosomeLength) / dataPanelWidth;

        // The test is true if we are zoomed in
        if (getReferenceFrame().getZoom() > 0) {

            double scale = getReferenceFrame().getScale();

            double origin = isDragging ? viewOrigin : getReferenceFrame().getOrigin();

            int start = (int) (origin / cytobandScale);
            double scaledDataPanelWidth = dataPanelWidth * scale;
            int span = (int) (scaledDataPanelWidth / cytobandScale);

            // Draw Cytoband current region viewer
            int height = (int) cytoRect.getHeight();
            g.setColor(Color.RED);
            int y = (int) (cytoRect.getY()) + CytobandRenderer.CYTOBAND_Y_OFFSET;
            currentRegionRect = new Rectangle(start - 2, y, span + 4, height);
            g.drawRect(start, y, span, height);
            g.drawRect(start - 1, (y - 1), span + 2, height + 2);
            g.drawRect(start - 2, (y - 2), span + 4, height + 4);
            if (span < 2) {
                g.drawRect(start - 2, (y - 2), span + 4, height + 4);
            }
        }
    }

    private void initMouseAdapter() {

        setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
        this.setToolTipText(
                "<html>Click anywhere on the chromosome<br/>to center view at that location.");


        MouseInputAdapter mouseAdapter = new MouseInputAdapter() {

            int lastMousePressX;

            public void mouseClicked(MouseEvent e) {
                if(currentCytobands == null) return;

                final int mouseX = e.getX();
                final int clickCount = e.getClickCount();
                WaitCursorManager.CursorToken token = WaitCursorManager.showWaitCursor();
                try {

                    double newLocation = cytobandScale * mouseX;
                    if (clickCount > 1) {
                        final int newZoom = getReferenceFrame().getZoom() + 1;
                        getReferenceFrame().zoomTo(newZoom, newLocation);
                    } else {
                        getReferenceFrame().centerOnLocation(newLocation);
                    }

                } finally {
                    WaitCursorManager.removeWaitCursor(token);
                }
            }

            @Override
            public void mousePressed(MouseEvent e) {
                lastMousePressX = e.getX();
            }

            @Override
            public void mouseReleased(MouseEvent e) {
                if(currentCytobands == null) return;

                if (isDragging) {
                    WaitCursorManager.CursorToken token = WaitCursorManager.showWaitCursor();
                    try {
                        getReferenceFrame().setOrigin(viewOrigin);
                        getReferenceFrame().recordHistory();
                    } finally {
                        WaitCursorManager.removeWaitCursor(token);
                    }
                }
                isDragging = false;
            }

            @Override
            public void mouseDragged(MouseEvent e) {

                if(currentCytobands == null) return;

                if (!isDragging && (currentRegionRect != null && currentRegionRect.contains(e.getPoint()))) {
                    isDragging = true;
                    viewOrigin = getReferenceFrame().getOrigin();
                }

                int w = getWidth();
                double scale = getReferenceFrame().getScale();

                int x = (int) Math.max(0, Math.min(e.getX(), w * (cytobandScale - scale)));
                int delta = x - lastMousePressX;
                if ((delta != 0) && (cytobandScale > 0)) {
                    viewOrigin = Math.max(0, Math.min(viewOrigin + delta * cytobandScale, w * (cytobandScale - scale)));
                    repaint();
                }
                lastMousePressX = x;

            }

            @Override
            public void mouseEntered(MouseEvent e) {
            }

            @Override
            public void mouseExited(MouseEvent e) {

            }
        };

        addMouseMotionListener(mouseAdapter);
        addMouseListener(mouseAdapter);
    }

    // TODO remove reference to IGV.theInstance

    private ReferenceFrame getReferenceFrame() {
        return frame;
    }
}
