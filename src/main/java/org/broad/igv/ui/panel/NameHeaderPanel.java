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
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.panel;

import org.broad.igv.Globals;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.IGVMouseInputAdapter;

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
