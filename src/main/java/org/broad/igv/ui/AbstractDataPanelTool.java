/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.ui;


import org.broad.igv.ui.panel.DataPanel;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.DataPanelTool;

import java.awt.*;
import java.awt.event.MouseEvent;

/**
 * @author eflakes
 */
abstract public class AbstractDataPanelTool implements DataPanelTool {

    private Cursor cursor;
    private String name;
    private DataPanel owner;

    public AbstractDataPanelTool(DataPanel owner, Cursor cursor) {
        this.owner = owner;
        this.cursor = cursor;
    }

    public Cursor getCursor() {
        return cursor;
    }

    public ReferenceFrame getReferenceFame() {
        return owner.getFrame();
    }

    public DataPanel getOwner() {
        return owner;
    }


    public void setName(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }

    public void mouseClicked(MouseEvent mouseEvent) {

    }

    public void mouseDragged(MouseEvent mouseEvent) {

    }

    public void mousePressed(MouseEvent mouseEvent) {

    }

    public void mouseReleased(MouseEvent mouseEvent) {
        
    }
}
