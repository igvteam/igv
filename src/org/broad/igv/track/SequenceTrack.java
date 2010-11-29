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

    /**
     * Field description
     */
    public static final int SEQUENCE_HEIGHT = 14;

    SequenceRenderer sequenceRenderer = new SequenceRenderer();
    ColorScale colorScale;
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

    public void setColorScale(ContinuousColorScale colorScale) {
        //To change body of implemented methods use File | Settings | File Templates.
    }


    public void setStatType(WindowFunction type) {

        // ignore
    }

    public WindowFunction getWindowFunction() {

        // ignore
        return null;
    }


    public void setRendererClass(Class rc) {

        // ignored
    }


    public Renderer getRenderer() {
        return sequenceRenderer;
    }


    @Override
    public int getHeight() {
        return SEQUENCE_HEIGHT * (showColorSpace ? 2 : 1);
    }

    public boolean isLogNormalized() {
        return true;
    }


    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {

        // TODO -- return sequence at this position
        return "";
    }


    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, ReferenceFrame frame) {
        return 0;
    }


    public boolean handleClick(int x, int y) {

        // Ignore
        return false;
    }


}
