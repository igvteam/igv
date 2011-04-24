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

package org.broad.igv.peaks;

import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.IGVMainFrame;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.*;
import java.io.IOException;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 * @date Apr 22, 2011
 */
public class PeakTrack extends AbstractTrack {

    enum ShadeOption {
        SCORE, FOLD_CHANGE
    }

    ;

    static PeakControlDialog controlDialog;
    private static float scoreThreshold = 0;
    private static ShadeOption shadeOption = ShadeOption.SCORE;

    int nTimePoints;
    Map<String, List<Peak>> peakMap = new HashMap();
    Renderer renderer = new PeakRenderer();
    private int bandHeight = 20;

    public PeakTrack(ResourceLocator locator, Genome genome) throws IOException {
        super(locator);
        height = bandHeight;
        loadPeaks(locator.getPath());
    }

    private void loadPeaks(String path) throws IOException {
        PeakParser parser = new PeakParser();
        List<Peak> peaks = parser.loadPeaks(path);
        nTimePoints = parser.getnTimePoints();
        TrackProperties props = parser.getTrackProperties();
        if (props != null) {
            setTrackProperties(props);
        }

        for (Peak peak : peaks) {
            String chr = peak.getChr();
            List<Peak> peakList = peakMap.get(chr);
            if (peakList == null) {
                peakList = new ArrayList();
                peakMap.put(chr, peakList);
            }
            peakList.add(peak);
        }

    }


    @Override
    public JPopupMenu getPopupMenu(TrackClickEvent te) {
        return new PeakTrackMenu(this);
    }

    public void render(RenderContext context, Rectangle rect) {

        List<Peak> peakList = peakMap.get(context.getChr());
        if (peakList == null) {
            return;
        }

        if (getDisplayMode() == Track.DisplayMode.EXPANDED) {
            Graphics2D borderGraphics = context.getGraphic2DForColor(Color.black);
            borderGraphics.drawLine(rect.x, rect.y, rect.x + rect.width, rect.y);
            borderGraphics.drawLine(rect.x, rect.y, rect.x + rect.width, rect.y);
            borderGraphics.drawLine(rect.x, rect.y + height, rect.x + rect.width, rect.y + height);
            borderGraphics.drawLine(rect.x, rect.y + height - 1, rect.x + rect.width, rect.y + height - 1);
            rect.y += 2;
            rect.height -= 4;
        }
        renderer.render(peakList, context, rect, this);
    }

    public Renderer getRenderer() {
        return renderer;
    }

    @Override
    public void setHeight(int height) {
        if (getDisplayMode() == Track.DisplayMode.COLLAPSED) {
            bandHeight = height;
        } else {
            bandHeight = (height - 6) / (nTimePoints + 1);
        }
        super.setHeight(height);
    }

    @Override
    public void setDisplayMode(DisplayMode mode) {
        super.setDisplayMode(mode);
        if (getDisplayMode() == Track.DisplayMode.COLLAPSED) {
            super.setHeight(bandHeight);
        } else if (getDisplayMode() == Track.DisplayMode.EXPANDED) {
            super.setHeight((nTimePoints + 1) * bandHeight + 6);
        }
    }


    // TODO -- the code below is an exact copy of code in DataTrack.   Refactor to share this.

    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {
        StringBuffer buf = new StringBuffer();
        Peak score = getLocusScoreAt(chr, position, frame);
        buf.append((score == null) ? "" : score.getValueString(position, getWindowFunction()));
        return buf.toString();
    }


    // TODO -- the code below is an exact copy of code in DataTrack.   Refactor to share this.

    private Peak getLocusScoreAt(String chr, double position, ReferenceFrame frame) {
        List<Peak> scores = peakMap.get(chr);

        if (scores == null) {
            return null;
        } else {
            // give a 2 pixel window, otherwise very narrow features will be missed.
            double bpPerPixel = frame.getScale();
            int buffer = (int) (2 * bpPerPixel);    /* * */
            return (Peak) FeatureUtils.getFeatureAt(position, buffer, scores);
        }
    }


    public static boolean controlDialogIsOpen() {
        return controlDialog != null && controlDialog.isVisible();
    }


    static synchronized void openControlDialog() {
        if (controlDialog == null) {
            controlDialog = new PeakControlDialog(IGV.getMainFrame());
        }
        controlDialog.setVisible(true);
    }


    public static float getScoreThreshold() {
        return scoreThreshold;
    }

    public static void setScoreThreshold(float t) {
        scoreThreshold = t;
    }

    public static ShadeOption getShadeOption() {
        return shadeOption;
    }

    public static void setShadeOption(ShadeOption shadeOption) {
        PeakTrack.shadeOption = shadeOption;
    }

}
