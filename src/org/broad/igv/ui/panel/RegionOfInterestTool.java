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

import org.broad.igv.feature.RegionOfInterest;
import org.broad.igv.ui.AbstractDataPanelTool;
import org.broad.igv.ui.IGV;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;

/**
 * @author eflakes
 */
public class RegionOfInterestTool extends AbstractDataPanelTool {

    Integer roiStart = null;
    JButton roiButton;

    public RegionOfInterestTool(DataPanel owner, JButton roiButton) {
        super(owner, Cursor.getPredefinedCursor(Cursor.CROSSHAIR_CURSOR));
        this.roiButton = roiButton;
        setName("Region of Interest");
    }

    public int getRoiStart() {
        return (roiStart == null ? 0 : roiStart.intValue());
    }

    /**
     * The mouse has been clicked.  Define one edge of the region of interest.
     */
    @Override
    public void mouseClicked(final MouseEvent e) {

        if (e.isPopupTrigger()) {
            return;
        }

        if(e.getClickCount() > 1) {
            return;
        }

        ReferenceFrame referenceFrame = this.getReferenceFame();

        if (e.getButton() == MouseEvent.BUTTON1 &&
                e.getClickCount() == 1) {

            Object chromosome = referenceFrame.getChromosome();

            // Allow Regions of Interest edits if ROI is enabled
            // and we have a valid Chromosome
            if (chromosome != null) {

                String chromosomeName = referenceFrame.getChrName();
                if (chromosomeName != null) {

                    int x = e.getX();

                    // Create a user Region of Interest
                    if (roiStart == null) {
                        roiStart = (int) referenceFrame.getChromosomePosition(x);
                        getOwner().repaint();
                    } else {

                        try {

                            int roiEnd = (int) referenceFrame.getChromosomePosition(x);
                            int start = Math.min(roiStart, roiEnd);
                            int end = Math.max(roiStart, roiEnd);

                            if (start == end) {
                                ++end;
                            }

                            // Create a Region of Interest
                            RegionOfInterest regionOfInterest =
                                    new RegionOfInterest(
                                            chromosomeName,
                                            start,
                                            end,
                                            null);
                            // TODO -- get this ugly reference out of here
                            IGV.getInstance().endROI();
                            IGV.getInstance().addRegionOfInterest(regionOfInterest);

                        } finally {
                            roiButton.setSelected(false);
                        }
                    }
                }

                IGV.getInstance().doRefresh();
            }
        }
    }

}
