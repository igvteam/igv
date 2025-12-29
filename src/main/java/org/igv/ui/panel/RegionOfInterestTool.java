/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.igv.ui.panel;

import org.igv.feature.RegionOfInterest;
import org.igv.ui.AbstractDataPanelTool;
import org.igv.ui.IGV;
import org.igv.ui.util.UIUtilities;

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

        if (e.getButton() == MouseEvent.BUTTON1 && e.getClickCount() == 1) {

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
                        UIUtilities.invokeOnEventThread(() -> getOwner().paintImmediately(getOwner().getBounds()));

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

                            IGV.getInstance().endROI();
                            IGV.getInstance().addRegionOfInterest(regionOfInterest);
                            IGV.getInstance().repaint();
                        } finally {
                            roiButton.setSelected(false);
                        }
                    }
                }

            }
        }
    }
}
