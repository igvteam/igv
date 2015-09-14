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

package org.broad.igv.cursor;


import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.Serializable;

/**
 * @author jrobinso
 *         Date: 1/26/14
 *         Time: 10:55 PM
 */
public class CursorTrackLabelPanel extends JComponent implements Serializable {

    final CursorModel model;
    final CursorTrack track;
    int sortDirection = 1;

    JButton sortButton;
    JRadioButton selectButton;
    CursorIdeogramPanel ideogramPanel;

    public CursorTrackLabelPanel(final CursorTrack track, final CursorModel model, final CursorMainPanel mainPanel) {

        this.track = track;
        this.model = model;
        this.sortButton = new JButton("Sort");
        this.sortButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                model.sortFrames(track, sortDirection);
                sortDirection *= -1;
                mainPanel.repaint();
            }
        });
        add(sortButton);

        selectButton = new JRadioButton();
        add(selectButton);
        mainPanel.addTrackSelectionButton(selectButton);

        selectButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(selectButton.isSelected()) {
                    model.setRegions(CursorUtils.createRegions(track));
                    mainPanel.repaint();

                }
            }
        });

//        ideogramPanel = new CursorIdeogramPanel();
//        ideogramPanel.addTrack(track);
//        ideogramPanel.setModel(model);
//        ideogramPanel.setDrawViewRect(false);
//        add(ideogramPanel);
    }


    @Override
    public void doLayout() {
        synchronized (this.getTreeLock()) {
            Rectangle bound = getBounds();
            sortButton.setBounds(0, 0, 50, bound.height);

            selectButton.setBounds(60, 0, bound.width - 40, bound.height);
        }
    }
}
