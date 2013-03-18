/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.panel;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.ui.AbstractDataPanelTool;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.event.DragStoppedEvent;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;

/**
 *
 */
public class PanTool extends AbstractDataPanelTool {

    private int previousYDirection = 0;    // Drag Directions: 1=up, 0=none 0r -1=down
    //private int lastMousePressedY;
    private int cumulativeDeltaX;
    private int cumulativeDeltaY;
    private Point lastMousePoint;
    //private JViewport viewport;
    private JScrollBar verticalScrollBar;
    private boolean isDragging = false;
    private Cursor dragCursor;

    private ReferenceFrame referenceFrame;

    public PanTool(DataPanel owner) {
        super(owner, Cursor.getDefaultCursor());
        this.dragCursor = IGV.fistCursor;
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
        return isDragging ? dragCursor : Cursor.getDefaultCursor();
    }

    public Point getLastMousePoint() {
        return lastMousePoint;
    }

    @Override
    public void mousePressed(final MouseEvent e) {

        if(e.isPopupTrigger()) {
            return;
        }
        
        lastMousePoint = e.getPoint();
        cumulativeDeltaX = 0;
        cumulativeDeltaY = 0;

    }

    public void setReferenceFrame(ReferenceFrame frame){
        this.referenceFrame = frame;
    }


    @Override
    public ReferenceFrame getReferenceFame() {
        if(referenceFrame != null) return referenceFrame;
        return super.getReferenceFame();
    }

    public void mouseReleased(final MouseEvent e) {

        if (isDragging) {
            isDragging = false;
            getReferenceFame().getEventBus().post(new DragStoppedEvent());
        }
        Component panel = (Component) e.getSource();
        panel.setCursor(getCursor());
    }


    @Override
    final public void mouseDragged(final MouseEvent e) {


        try {
            Component panel = (Component) e.getSource();            
            panel.setCursor(dragCursor);
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
                    getReferenceFame().shiftOriginPixels(deltaX);
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

}
