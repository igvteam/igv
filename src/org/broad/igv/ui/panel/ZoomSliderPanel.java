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
 * ZoomSliderPanel.java
 *
 * Created on September 25, 2007, 12:03 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package org.broad.igv.ui.panel;

import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.event.ViewChange;
import org.broad.igv.ui.util.IconFactory;

import javax.swing.*;
import javax.swing.event.MouseInputAdapter;
import java.awt.*;
import java.awt.event.MouseEvent;

/**
 * @author jrobinso
 */
public class ZoomSliderPanel extends JPanel {

    static Color TICK_GRAY = new Color(90, 90, 90);
    static Color TICK_BLUE = new Color(25, 50, 200);

    //double imageScaleFactor = 0.8;
    Image slider;
    Image zoomPlus;
    Image zoomMinus;
    Rectangle zoomPlusRect;
    Rectangle zoomMinusRect;
    Rectangle[] zoomLevelRects;
    /**
     * Should correspond to "maxZoomLevel" in class referenceFrame.
     */
    int numZoomLevels = 25;


    private int minZoomLevel = 0;

    /**
     * Set the allowed zoom level, user cannot zoom out past this level
     *
     * @param minZoomLevel
     */
    public void setMinZoomLevel(int minZoomLevel){
        this.minZoomLevel = minZoomLevel;
    }

    private static final Color TRANSPARENT_GRAY = new Color(200, 200, 200, 150);
    private ReferenceFrame referenceFrame;

    public ZoomSliderPanel(){
        this(null);
    }

    /**
     * @param referenceFrame The ReferenceFrame whose zoom level this panel will control
     */
    public ZoomSliderPanel(ReferenceFrame referenceFrame) {
        this.referenceFrame = referenceFrame;
        slider = IconFactory.getInstance().getIcon(IconFactory.IconID.SLIDER).getImage();
        zoomPlus = IconFactory.getInstance().getIcon(IconFactory.IconID.ZOOM_PLUS).getImage();
        zoomMinus = IconFactory.getInstance().getIcon(IconFactory.IconID.ZOOM_MINUS).getImage();
        zoomLevelRects = new Rectangle[numZoomLevels];

        setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));

        init();
    }

    private void updateTickCount() {
        int tmp = getViewContext().getMaxZoom() + 1;
        if (tmp != numZoomLevels) {
            numZoomLevels = tmp;
            zoomLevelRects = new Rectangle[numZoomLevels];
        }

    }


    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        updateTickCount();
        //if (this.isEnabled()) {
        paintHorizontal(g);
        //}

    }

    protected void paintHorizontal(Graphics g) {

        Graphics2D transGraphics = (Graphics2D) g.create();
        transGraphics.setColor(TRANSPARENT_GRAY);

        int buttonWidth = zoomPlus.getWidth(null);
        int buttonHeight = zoomPlus.getHeight(null);

        Insets insets = getInsets();
        int panelWidth = getWidth() - insets.left - insets.right;
        int panelHeight = getHeight() - insets.top - insets.bottom;


        final boolean enabled = isEnabled();
        g.setColor(enabled ? Color.BLACK : Color.LIGHT_GRAY);
        double x = insets.left;

        double xStep = ((double) (panelWidth - 2 * buttonWidth - 10)) / (numZoomLevels);

        int y = insets.top + (panelHeight - buttonHeight) / 2;
        g.drawImage(zoomMinus, (int) x, y, null);
        zoomMinusRect = new Rectangle((int) x, y, buttonWidth, buttonHeight);

        if (!isEnabled()) {
            transGraphics.fill(zoomMinusRect);
        }

        x += 5 + buttonWidth;

        int lastX = (int) (x - xStep);
        for (int i = 0; i < numZoomLevels; i++) {
            Rectangle zoomRect = new Rectangle((int) x, y, (int) (x - lastX), buttonHeight);
            int xLine = (int) (x + xStep / 2);
            g.drawLine(xLine, y + 3, xLine, y + buttonHeight - 4);
            zoomLevelRects[i] = zoomRect;
            lastX = (int) x;
            x += xStep;
        }

        x += 5;

        y = insets.top + panelHeight / 2 - 1;
        //g.drawLine(xTop, y, xBottom, y);

        y = insets.top + (panelHeight - buttonHeight) / 2;
        //if (isEnabled()) {
        g.drawImage(zoomPlus, (int) x, y, null);
        //}
        zoomPlusRect = new Rectangle((int) x, y, buttonWidth, buttonWidth);

        if (!isEnabled()) {
            transGraphics.fill(zoomPlusRect);
        }

        // Draw current level -- zoomIndex is the zoom level + 1. 

        int zoom = (toolZoom >= 0 ? toolZoom : getViewContext().getAdjustedZoom());

        if (enabled) {
            if (zoom >= 0 && zoom < zoomLevelRects.length) {
                Rectangle rect = zoomLevelRects[zoom];

                g.setColor(TICK_BLUE);
                g.fill3DRect(
                        (int) (rect.getX() + rect.getWidth() / 2) - 3,
                        (int) rect.getY(),
                        6,
                        (int) rect.getHeight(),
                        true);

                //y =  (int) (rect.getY() + (rect.getHeight() - slider.getHeight(null)) / 2);
                // temporary hack
                //if(zoomIndex == 12) y += 15;
                //g.drawImage(slider, x + 1, y, null);
            }
        }
        transGraphics.dispose();
    }

    int setZoom(MouseEvent e) {

        if (zoomPlusRect.contains(e.getX(), e.getY())) {
            toolZoom++;
        } else if (zoomMinusRect.contains(e.getX(), e.getY()) && toolZoom > minZoomLevel) {
            toolZoom--;
        } else {
            for (int i = 0; i < zoomLevelRects.length; i++) {
                Rectangle rect = zoomLevelRects[i];
                if (rect.contains(e.getX(), e.getY()) && i >= minZoomLevel) {
                    toolZoom = i;
                }
            }
        }
        return toolZoom;
    }


    private ReferenceFrame getViewContext() {
        if(referenceFrame == null) return FrameManager.getDefaultFrame();
        return referenceFrame;
    }

    int toolZoom = -1;

    private void init() {

        MouseInputAdapter mouseAdapter = new MouseInputAdapter() {

            int lastMousePressX = 0;

            @Override
            public void mouseExited(MouseEvent e) {

            }

            @Override
            public void mouseClicked(MouseEvent e) {

            }

            @Override
            public void mouseMoved(MouseEvent e) {
            }

            @Override
            public void mousePressed(MouseEvent e) {
                if (!isEnabled()) {
                    return;
                }
                toolZoom = Math.max(0, getViewContext().getAdjustedZoom());
            }

            @Override
            public void mouseReleased(MouseEvent e) {
                if (!isEnabled()) {
                    return;
                }
                //Sometimes the zoom doesn't change, don't need to do anything in that case
                int oldToolZoom = toolZoom;
                int diff = setZoom(e) - oldToolZoom;
                if(diff == 0) {
                    toolZoom = -1;
                    return;
                }

                repaint();

                int effectiveZoom = toolZoom + getViewContext().getMinZoom();

                ViewChange.ZoomCause event = new ViewChange.ZoomCause(effectiveZoom);
                getViewContext().getEventBus().post(event);
                toolZoom = -1;

            }


            @Override
            public void mouseDragged(MouseEvent e) {
                // Dragging zoom tool is disable.  Generates too many
                // repaint events.
                setZoom(e);
                repaint();
            }
        };

        addMouseMotionListener(mouseAdapter);
        addMouseListener(mouseAdapter);
    }
}
