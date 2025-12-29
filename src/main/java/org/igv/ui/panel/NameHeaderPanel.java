/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.igv.ui.panel;

import org.igv.Globals;
import org.igv.ui.IGV;
import org.igv.ui.util.IGVMouseInputAdapter;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;

/**
 * Name header panel.  Displayed in upper left corner
 *
 * @author jrobinso
 */
public class NameHeaderPanel extends JPanel implements Paintable {


    private final boolean darkMode;

    public NameHeaderPanel() {

        this.darkMode = Globals.isDarkMode();
        
        // Clicking on panel will clear all track selections
        this.addMouseListener(new IGVMouseInputAdapter() {
            @Override
            public void igvMouseClicked(MouseEvent mouseEvent) {
                IGV.getInstance().clearSelections();
                IGV.getInstance().repaintNamePanels();
            }
        });
    }

    public void paintOffscreen(Graphics2D g, Rectangle rect, boolean batch) {
        paintComponent(g);
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);

        if(darkMode){
            setBackground(UIManager.getColor("Panel.background"));
        }

    }

    @Override
    public int getSnapshotHeight(boolean batch) {
        return getHeight();
    }
}
