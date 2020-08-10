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
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.event.IGVEventObserver;
import org.broad.igv.event.ViewChange;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.Cytoband;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.renderer.CytobandRenderer;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.ui.util.IGVMouseInputAdapter;

import javax.swing.*;
import javax.swing.event.MouseInputAdapter;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.util.List;

/**
 * @author jrobinso
 */
public class CytobandPanel extends JPanel implements IGVEventObserver {

    private static int fontHeight = 10;
    private static int bandHeight = 10;
    private static String fontFamilyName = "Lucida Sans";
    private boolean isDragging = false;

    /**
     * left genomic coordinate of region in view (base pairs)
     */
    private double viewOrigin;

    /**
     * right genomic coordinate of region in view (base pairs)
     */

    private double viewEnd;

    /**
     * Scale in base-pairs per pixel == chromosome length / panel width
     */
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
        this.viewOrigin = frame.getOrigin();
        this.viewEnd = frame.getEnd();

        FontManager.getFont(fontHeight);
        setFont(new Font(fontFamilyName, Font.BOLD, fontHeight));
        if (mouseable) {
            initMouseAdapter();
        }
        cytobandRenderer = (new CytobandRenderer());
        IGVEventBus.getInstance().subscribe(ViewChange.class, this);

    }


    @Override
    protected void paintComponent(Graphics g) {

        super.paintComponent(g);

        if (PreferencesManager.getPreferences().getAntiAliasing()) {
            ((Graphics2D) g).setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        }

        if (frame.getChrName().equals(Globals.CHR_ALL) || getWidth() < 10) {
            return;
        }

        int dataPanelWidth = getWidth();
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

        int chromosomeLength = getReferenceFrame().getMaxCoordinate();
        cytobandScale = ((double) chromosomeLength) / dataPanelWidth;

        // The test is true if we are zoomed in
        if (getReferenceFrame().getZoom() > 0) {

            double origin = viewOrigin;
            double end = viewEnd;

            int pixelStart = (int) (origin / cytobandScale);
            int pixelEnd = (int) (end / cytobandScale);
            int pixelSpan = Math.max(0, pixelEnd - pixelStart);

            // Draw Cytoband current region viewer
            int height = (int) cytoRect.getHeight();
            g.setColor(Color.RED);
            int y = (int) (cytoRect.getY()) + CytobandRenderer.CYTOBAND_Y_OFFSET;
            currentRegionRect = new Rectangle(pixelStart - 2, y, pixelSpan + 4, height);
            g.drawRect(pixelStart, y, pixelSpan, height);
            g.drawRect(pixelStart - 1, (y - 1), pixelSpan + 2, height + 2);
            g.drawRect(pixelStart - 2, (y - 2), pixelSpan + 4, height + 4);
            if (pixelSpan < 2) {
                g.drawRect(pixelStart - 2, (y - 2), pixelSpan + 4, height + 4);
            }
        }
    }

    private void initMouseAdapter() {

        setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
        this.setToolTipText(
                "<html>Click anywhere on the chromosome<br/>to center view at that location.");


        MouseInputAdapter mouseAdapter = new IGVMouseInputAdapter() {

            int lastMousePressX;

            public void igvMouseClicked(MouseEvent e) {
                if (currentCytobands == null) return;

                final int mouseX = e.getX();
                final int clickCount = e.getClickCount();
                WaitCursorManager.CursorToken token = WaitCursorManager.showWaitCursor();
                try {

                    double newLocation = cytobandScale * mouseX;
                    if (clickCount > 1) {
                        final int newZoom = getReferenceFrame().getZoom() + 1;
                        getReferenceFrame().doSetZoomCenter(newZoom, newLocation);
                    } else {
                        getReferenceFrame().centerOnLocation(newLocation);
                    }

                    ViewChange result = ViewChange.Result();
                    result.setRecordHistory(true);
                    IGVEventBus.getInstance().post(result);

                } finally {
                    WaitCursorManager.removeWaitCursor(token);
                }
            }

            @Override
            public void mousePressed(MouseEvent e) {
                super.mousePressed(e);
                lastMousePressX = e.getX();
            }

            @Override
            public void mouseReleased(MouseEvent e) {
                super.mouseReleased(e);
                if (currentCytobands == null) return;

                if (isDragging) {
                    WaitCursorManager.CursorToken token = WaitCursorManager.showWaitCursor();
                    try {
                        getReferenceFrame().setOrigin(viewOrigin);
                        getReferenceFrame().recordHistory();
                    } finally {
                        WaitCursorManager.removeWaitCursor(token);
                    }
                    ViewChange result = ViewChange.Result();
                    result.setRecordHistory(true);
                    IGVEventBus.getInstance().post(result);
                }
                isDragging = false;
            }

            @Override
            public void mouseDragged(MouseEvent e) {
                if (currentCytobands == null) return;
                isDragging = true;
                int x = e.getX();
                int delta = x - lastMousePressX;

                if ((delta != 0) && (cytobandScale > 0)) {
                    // Constrain to bounds of chromosome
                    double chrLength = CytobandPanel.this.frame.getChromosomeLength();
                    double deltaBP = Math.min( Math.max(-viewOrigin, delta * cytobandScale), chrLength - viewEnd);
                    viewOrigin += deltaBP;
                    viewEnd += deltaBP;
                    // TODO Constrain to chromosome bounds?
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

    private ReferenceFrame getReferenceFrame() {
        return frame;
    }

    @Override
    public void receiveEvent(Object e) {
        if (e instanceof ViewChange) {
            ViewChange event = (ViewChange) e;
            if (event.type == ViewChange.Type.ChromosomeChange || event.type == ViewChange.Type.LocusChange) {
                this.viewOrigin = frame.getOrigin();
                this.viewEnd = frame.getEnd();
            }
        }
    }
}
