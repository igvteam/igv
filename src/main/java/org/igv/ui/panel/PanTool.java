/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.igv.ui.panel;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.commons.math3.stat.StatUtils;
import org.igv.ui.AbstractDataPanelTool;
import org.igv.ui.IGV;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Map;

/**
 *
 */
public class PanTool extends AbstractDataPanelTool {

    private static int throttleTimeMS = 50;

    private static double [] repaintTimes = new double[100];
    private static int repaintTimeCount = 0;

    private int previousYDirection = 0;    // Drag Directions: 1=up, 0=none 0r -1=down
    //private int lastMousePressedY;
    private int cumulativeDeltaX;
    private int cumulativeDeltaY;
    private Point lastMousePoint;
    //private JViewport viewport;
    private JScrollBar verticalScrollBar;
    private boolean isDragging = false;
    private long lastDragEventTime = 0;

    private ReferenceFrame referenceFrame;

    public PanTool(DataPanel owner) {
        super(owner, Cursor.getDefaultCursor());
        setName("Pan");
        if (owner != null) {
            verticalScrollBar = owner.getVerticalScrollbar();
            Container parentContainer = owner.getParent();
            if (parentContainer != null) {
                Container parentOfParent = parentContainer.getParent();
                if ((parentOfParent != null) && (parentOfParent instanceof JViewport)) {
                    //viewport = (JViewport) parentOfParent;
                }
            }
        }
    }


    @Override
    public Cursor getCursor() {
        return Cursor.getDefaultCursor();
    }

    public Point getLastMousePoint() {
        return lastMousePoint;
    }

    @Override
    public void mousePressed(final MouseEvent e) {

        if (e.isPopupTrigger()) {
            return;
        }

        lastMousePoint = e.getPoint();
        cumulativeDeltaX = 0;
        cumulativeDeltaY = 0;

    }

    public void setReferenceFrame(ReferenceFrame frame) {
        this.referenceFrame = frame;
    }


    @Override
    public ReferenceFrame getReferenceFame() {
        if (referenceFrame != null) {
            return referenceFrame;
        }
        return super.getReferenceFame();
    }

    public void mouseReleased(final MouseEvent e) {

        if (isDragging) {
            isDragging = false;
            lastDragEventTime = 0;
            getReferenceFame().dragStopped();
        }
        Component panel = (Component) e.getSource();
        panel.setCursor(getCursor());
    }


    @Override
    final public void mouseDragged(final MouseEvent e) {

        long currentTime = System.currentTimeMillis();

        if ((currentTime - lastDragEventTime) < throttleTimeMS) {
            return;
        }

        //System.out.println("T=" + ctr);
        lastDragEventTime = currentTime;

        try {
            Component panel = (Component) e.getSource();
            if (lastMousePoint == null) {
                lastMousePoint = e.getPoint();
                return;
            }

            if (!isDragging && e.getPoint().distance(lastMousePoint) < 2) {
                return;
            } else {
                isDragging = true;

                int deltaX = lastMousePoint.x - e.getX();
                int deltaY = lastMousePoint.y - e.getY();
                cumulativeDeltaX += Math.abs(deltaX);
                cumulativeDeltaY += Math.abs(deltaY);


                // Test for horizontal vs vertical panning.
                if (cumulativeDeltaX > cumulativeDeltaY) {

                    // Horizontal scrolling
                    getReferenceFame().shiftOriginPixelsPanning(deltaX);
                } else {
                    // Vertical Scrolling 
                    int totalYChange = (int) (lastMousePoint.getY() - e.getY());
                    if (totalYChange != 0) {

                        // This section handles false drag direction changes
                        int currentYDirection = 0;

                        // Figure out the current drag direction
                        currentYDirection = totalYChange / Math.abs(totalYChange);


                        // If the previous direction is 0 we were not moving before
                        if (previousYDirection != 0) {

                            // See if we changed direction
                            boolean changedYDirection = currentYDirection != previousYDirection;
                            if (!changedYDirection) {

                                // Don't save lastMousePressedPoint because may
                                // be incorrect (this is the problem we are
                                // solving with the direction flag) instead
                                // we'll just check the next drag Point to be
                                // sure of the correct direction.
                                previousYDirection = currentYDirection;

                                // If we have a vertical scrollbar use it to move
                                if (verticalScrollBar != null) {
                                    int adjustedScrollbarValue = verticalScrollBar.getValue();
                                    adjustedScrollbarValue += totalYChange;
                                    verticalScrollBar.setValue(adjustedScrollbarValue);
                                }
                            }
                        }
                        previousYDirection = currentYDirection;

                    }
                }
            }
        } finally {
            lastMousePoint = e.getPoint();    // Always save the last Point
        }
    }

    /**
     * Adaptively throttle the drag time based on repaint times.  Its pointless to drag faster than we can repaint.
     *
     * @param t
     */
    public static void repaintTime(double t) {
            repaintTimes[repaintTimeCount++] = t;
            if (repaintTimeCount == repaintTimes.length) {
                double p = StatUtils.percentile(repaintTimes, 95);
                throttleTimeMS = Math.max(5, Math.min(200, (int) p));
                repaintTimeCount = 0;
            }

    }

}
