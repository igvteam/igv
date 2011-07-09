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


package org.broad.igv.track;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.Strand;
import org.broad.igv.renderer.*;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.IGVPopupMenu;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;


/**
 * @author jrobinso
 */
public class SequenceTrack extends AbstractTrack {

    private static Logger log = Logger.getLogger(SequenceTrack.class);

    private static final int SEQUENCE_HEIGHT = 14;

    private static String NAME = "Sequence";

    private SequenceRenderer sequenceRenderer = new SequenceRenderer();

    //should translated aminoacids be shown below the sequence?
    private boolean shouldShowTranslation = true;

    //is sequence visible (zoomed in far enough, etc)
    private boolean sequenceVisible = false;

    Strand strand = Strand.POSITIVE;

    /**
     * If true show sequence in "color space"  (for SOLID alignments).  Currently not implemented.
     */
    private boolean showColorSpace = false;
    private Rectangle arrowRect;

    public SequenceTrack(String name) {
        super(name);
        setSortable(false);
        shouldShowTranslation = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_SEQUENCE_TRANSLATION);

    }


    @Override
    public void renderName(Graphics2D graphics, Rectangle trackRectangle, Rectangle visibleRectangle) {
        Font font= FontManager.getFont(fontSize);
        if (sequenceVisible) {
            graphics.setFont(font);
            graphics.drawString(NAME, trackRectangle.x + 5, trackRectangle.y + 12);

            int rx = trackRectangle.x + trackRectangle.width - 20;
            arrowRect = new Rectangle(rx, trackRectangle.y + 2, 15, 10);
            drawArrow(graphics);
        }
    }

    private void drawArrow(Graphics2D graphics) {
        if (strand == Strand.POSITIVE) {
            GraphicUtils.drawHorizontalArrow(graphics, arrowRect, true);
        } else {
            GraphicUtils.drawHorizontalArrow(graphics, arrowRect, false);
        }
    }

    /**
     * Render the sequence, and optionally the 3 frame translation table
     *
     * @param context
     * @param rect
     */
    public void render(RenderContext context, Rectangle rect) {
        // Are we zoomed in far enough to show the sequence?  Scale is
        // in BP / pixel,  need at least 1 pixel  per bp in order to show sequence.

        int resolutionThreshold =  PreferenceManager.getInstance().getAsInt(PreferenceManager.MAX_SEQUENCE_RESOLUTION);
        // TODO -- this should be calculated from a "rescale" event
        boolean visible = FrameManager.getMinimumScale() < resolutionThreshold &&
                !context.getChr().equals(Globals.CHR_ALL);

        if (visible != sequenceVisible) {
            sequenceVisible = visible;
            IGV.getInstance().doRefresh();
        }
        if (sequenceVisible) {
            sequenceRenderer.setStrand(strand);
            sequenceRenderer.draw(context, rect, showColorSpace, shouldShowTranslation, resolutionThreshold);
        }
    }


    @Override
    public int getHeight() {
        return sequenceVisible ? SEQUENCE_HEIGHT + (showColorSpace ? SEQUENCE_HEIGHT : 0) +
                (shouldShowTranslation ? SequenceRenderer.TranslatedSequenceDrawer.TOTAL_HEIGHT : 0) :
                0;
    }

    @Override
    public int getPreferredHeight() {
        return getHeight();
    }


    @Override
    public boolean handleDataClick(TrackClickEvent e) {

        MouseEvent evt = e.getMouseEvent();
        setShouldShowTranslation(!shouldShowTranslation);
        Object source = e.getMouseEvent().getSource();
        if (source instanceof JComponent) {
            repaint();

        }
        return true;
    }

    @Override
    public void handleNameClick(final MouseEvent e) {
        if (arrowRect != null && arrowRect.contains(e.getPoint())) {
            flipStrand(e.getSource());
        }

    }

    private void flipStrand(Object source) {
        strand = (strand == Strand.POSITIVE ? Strand.NEGATIVE : Strand.POSITIVE);
        repaint();
        IGV.getInstance().getTrackManager().clearSelections();
    }

    public void setShouldShowTranslation(boolean shouldShowTranslation) {
        this.shouldShowTranslation = shouldShowTranslation;
        // Remember this choice
        PreferenceManager.getInstance().put(PreferenceManager.SHOW_SEQUENCE_TRANSLATION, shouldShowTranslation);
    }


    /**
     * Override to return a specialized popup menu
     *
     * @return
     */
    @Override
    public IGVPopupMenu getPopupMenu(final TrackClickEvent te) {

        IGVPopupMenu menu = new IGVPopupMenu();

        JMenuItem m1 = new JMenuItem("Flip strand");
        m1.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                flipStrand(te.getMouseEvent().getSource());

            }
        });

        final JCheckBoxMenuItem m2 = new JCheckBoxMenuItem("Show translation");
        m2.setSelected(shouldShowTranslation);
        m2.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                setShouldShowTranslation(m2.isSelected());
                repaint();
                IGV.getInstance().getTrackManager().clearSelections();
            }
        });

        menu.add(m1);
        menu.add(m2);

        return menu;
    }


    private void repaint() {
        // TODO -- what's really needed is a repaint of all panels the sequence track intersects
        IGV.getMainFrame().repaint();
    }

    // SequenceTrack does not expose its renderer

    public Renderer getRenderer() {
        return null;
    }


}
