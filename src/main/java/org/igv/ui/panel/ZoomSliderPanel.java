/*
 * ZoomSliderPanel.java
 *
 * Created on September 25, 2007, 12:03 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package org.igv.ui.panel;

import org.igv.Globals;
import org.igv.track.Track;
import org.igv.ui.IGV;
import org.igv.ui.util.IGVMouseInputAdapter;
import org.igv.ui.util.IconFactory;

import javax.swing.*;
import javax.swing.event.MouseInputAdapter;

import java.util.List;
import java.awt.*;
import java.awt.event.MouseEvent;

/**
 * @author jrobinso
 */
public class ZoomSliderPanel extends JPanel {

    private static final Color TRANSPARENT_GRAY = new Color(200, 200, 200, 150);
    private static final Color TRANSPARENT_BLUE = new Color(27, 96, 246, 25);
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
    private ReferenceFrame referenceFrame;
    private boolean darkMode;
    private Color markColor;
    private Color tickColor;


    /**
     * Set the allowed zoom level, user cannot zoom out past this level
     *
     * @param minZoomLevel
     */
    public void setMinZoomLevel(int minZoomLevel){
        this.minZoomLevel = minZoomLevel;
    }



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
        int tmp = Math.max(0, getReferenceFrame().getMaxZoom() + 1);
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
        g.setColor(enabled ? this.tickColor : Color.LIGHT_GRAY);
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

        int zoom = getReferenceFrame().getAdjustedZoom();

        if (enabled) {
            if (zoom >= 0 && zoom < zoomLevelRects.length) {
                Rectangle rect = zoomLevelRects[zoom];

                g.setColor(this.markColor);
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

        paintVisibilityThresholds(transGraphics);
        transGraphics.dispose();
    }

    //Adds an indicator of which zoom level will display reads/features
    private void paintVisibilityThresholds(final Graphics2D transGraphics) {
        if(numZoomLevels > 1) {
            List<Integer> visibilityThresholds = IGV.getInstance().getAllTracks().stream()
                    .map(Track::getVisibilityWindow)
                    .filter(i -> i > 0)
                    .sorted()
                    .distinct()
                    .map(threshold -> this.getReferenceFrame().calculateZoom(0, threshold))
                    .filter(z -> z > 1)
                    .toList();

            transGraphics.setColor(TRANSPARENT_BLUE);
            Rectangle maxZoom = zoomLevelRects[zoomLevelRects.length - 1];
            for (Integer window : visibilityThresholds) {
                final Rectangle currentLevel = zoomLevelRects[window];
                Rectangle windowBox = new Rectangle(currentLevel.x, currentLevel.y,
                        maxZoom.x + maxZoom.width - currentLevel.x, currentLevel.height);
                transGraphics.fill(windowBox);
            }
        }
    }

    void setZoom(MouseEvent e) {

        if (zoomPlusRect.contains(e.getX(), e.getY())) {
            getReferenceFrame().doZoomIncrement(1);
        } else if (zoomMinusRect.contains(e.getX(), e.getY())) {
            getReferenceFrame().doZoomIncrement(-1);
        } else {
            for (int i = 0; i < zoomLevelRects.length; i++) {
                Rectangle rect = zoomLevelRects[i];
                if (rect.contains(e.getX(), e.getY()) && i >= minZoomLevel) {
                    getReferenceFrame().setAdjustedZoom(i);
                    break;
                }
            }
        }
    }


    private ReferenceFrame getReferenceFrame() {
        if(referenceFrame == null) return FrameManager.getDefaultFrame();
        return referenceFrame;
    }


    private void init() {

        this.darkMode = Globals.isDarkMode();
        this.markColor = darkMode ? Color.CYAN : TICK_BLUE;
        this.tickColor = darkMode ? Color.white : Color.black;

        MouseInputAdapter mouseAdapter = new IGVMouseInputAdapter() {

            int lastMousePressX = 0;

            @Override
            public void mousePressed(MouseEvent e) {
                super.mousePressed(e);
                if (!isEnabled()) {
                    return;
                }
                //toolZoom = Math.max(0, getReferenceFrame().getAdjustedZoom());
            }

            @Override
            public void mouseReleased(MouseEvent e) {
                super.mouseReleased(e);
                if (!isEnabled()) {
                    return;
                }
                setZoom(e);
                repaint();
            }


            @Override
            public void mouseDragged(MouseEvent e) {
                // Dragging zoom tool is disable.  Generates too many
                // repaint events.
                //setZoom(e);
                //repaint();
            }
        };

        addMouseMotionListener(mouseAdapter);
        addMouseListener(mouseAdapter);
    }
}
