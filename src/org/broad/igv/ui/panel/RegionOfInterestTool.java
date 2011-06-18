/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
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
