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
