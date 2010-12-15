/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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


package org.broad.igv.track;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.Globals;
import org.broad.igv.renderer.ColorScale;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.renderer.SequenceRenderer;
import org.broad.igv.ui.panel.ReferenceFrame;

import java.awt.*;


/**
 * @author jrobinso
 */
public class SequenceTrack extends AbstractTrack {


    private static final int SEQUENCE_HEIGHT = 14;

    private SequenceRenderer sequenceRenderer = new SequenceRenderer();

    /**
     * If true show sequence in "color space"  (for SOLID alignments).  Currently not implemented, should always be
     * false.
     */
    private boolean showColorSpace = false;

    public SequenceTrack(String name) {
        super(name);
    }

    /**
     * Method description
     *
     * @param context
     * @param rect
     */
    public void render(RenderContext context, Rectangle rect) {

        // Are we zoomed in far enough to show the sequence?  Scale is
        // in BP / pixel,  need at least 1 pixel.
        if (context.getScale() < 1 && !context.getChr().equals(Globals.CHR_ALL)) {
            sequenceRenderer.draw(context, rect, showColorSpace);
        }

    }


    @Override
    public int getHeight() {
        return SEQUENCE_HEIGHT * (showColorSpace ? 2 : 1);
    }


    //----------------------------------------------------------------------------
    // Methods belowo are required for the Track interface, but aren't
    // meaningful here.  Obviously some refactoring is in order to reduce
    // the number of required methods.

    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {
        return null;
    }


    public void setColorScale(ContinuousColorScale colorScale) {
        // Required method for track interface, ignore
    }


    public void setStatType(WindowFunction type) {
        // Required method for track interface, ignore
    }

    public WindowFunction getWindowFunction() {
        // Required method for track interface, ignore
        return null;
    }


    public void setRendererClass(Class rc) {
        // Required method for track interface, ignore
    }


    public Renderer getRenderer() {
        // Required method for track interface, ignore
        return null;
    }


    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, ReferenceFrame frame) {
        // Required method for track interface, ignore
        return 0;
    }


    public boolean handleClick(int x, int y) {
        // Ignore
        return false;
    }


    public boolean isLogNormalized() {
        // Required method for track interface, ignore
        return true;
    }


}
