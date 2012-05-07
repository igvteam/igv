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
package org.broad.igv.sam;

import org.broad.igv.track.RenderContext;

import java.awt.*;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 */
public class BedRenderer implements FeatureRenderer {

    public void renderAlignments(List<Alignment> alignments, RenderContext context,
                                 Rectangle rowRect, Rectangle inputRect, AlignmentTrack.RenderOptions renderOptions, boolean leaveMargin,
                                 Map<String, Color> selectedReadNames) {

        double origin = context.getOrigin();
        double locScale = context.getScale();

        if ((alignments != null) && (alignments.size() > 0)) {


            for (Alignment alignment : alignments) {

                // Compute the start and dend of the alignment in pixels
                double pixelStart = ((alignment.getStart() - origin) / locScale);
                double pixelEnd = ((alignment.getEnd() - origin) / locScale);

                // If the any part of the feature fits in the track rectangle draw  it
                if ((pixelEnd >= rowRect.getX()) && (pixelStart <= rowRect.getMaxX())) {

                    Graphics2D g = context.getGraphic2DForColor(Color.GRAY);

                    // If the alignment is 3 pixels or less,  draw alignment as a single block,
                    // further detail would not be seen and just add to drawing overhead
                    int w = Math.max(1, (int) (pixelEnd - pixelStart));
                    int h = (int) Math.max(1, rowRect.getHeight() - 2);
                    int y = (int) (rowRect.getY() + (rowRect.getHeight() - h) / 2);
                    g.fillRect((int) pixelStart, y, w, h);


                }
            }

        }

    }
}
