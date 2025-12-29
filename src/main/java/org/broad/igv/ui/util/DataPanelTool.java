/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.ui.util;

import java.awt.*;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

/**
 * Defines an interface for a "tool" to be used in the data panel.  A tool handles a subset of mouse events for
 * the panel.  
 */
public interface DataPanelTool {

    String getName();

    Cursor getCursor();

    void mousePressed(MouseEvent mouseEvent);

    void mouseReleased(MouseEvent mouseEvent);

    void mouseClicked(MouseEvent mouseEvent);

    void mouseDragged(MouseEvent mouseEvent);

}
