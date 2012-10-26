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
package org.broad.igv.maf;

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.maf.MAFTile.MASequence;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.track.AbstractTrack;
import org.broad.igv.track.RegionScoreType;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.TrackClickEvent;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;

/**
 * @author jrobinso
 */
public class MultipleAlignmentTrack extends AbstractTrack {

    public static final int margin = 5;
    private static Logger log = Logger.getLogger(MultipleAlignmentTrack.class);
    private static int EXPANDED_HEIGHT = 14;
    private static int GAPS_HEIGHT = 25;
    //List<Rectangle> featureRects = new ArrayList();
    MAFManager mgr;
    MAFRenderer renderer = new MAFRenderer();
    // A hack until full MAF track is implemented.
    Rectangle visibleNameRect;

    Genome genome;

    /**
     * Map of chr alias -> "true" chromosome names.
     */
    private HashMap<String, String> chrMappings;


    public MultipleAlignmentTrack(ResourceLocator locator, Genome genome) throws IOException {
        super(locator);
        this.genome = genome;
        this.mgr = new MAFManager(locator, genome);

        List<String> mafChrNames = mgr.getChrNames();
        if (mafChrNames != null) {
            chrMappings = new HashMap();
            for (String mafChr : mafChrNames) {
                String chr = genome.getChromosomeAlias(mafChr);
                chrMappings.put(chr, mafChr);
            }
        }
    }

    @Override
    public int getHeight() {
        return GAPS_HEIGHT + (mgr.getSelectedSpecies().size() + 1) * EXPANDED_HEIGHT;
    }

    @Override
    public void renderName(Graphics2D g2D, Rectangle trackRectangle, Rectangle visibleRectangle) {

        this.visibleNameRect = trackRectangle;
        if (isSelected()) {
            g2D.setBackground(Color.LIGHT_GRAY);
        } else {
            g2D.setBackground(Color.WHITE);
        }

        Rectangle rect = new Rectangle(trackRectangle);
        g2D.clearRect(rect.x, rect.y, rect.width, rect.height);

        Font font = FontManager.getFont(fontSize);
        g2D.setFont(font);

        int y = trackRectangle.y;

        rect.height = GAPS_HEIGHT;
        rect.y = y;
        //text, rect.x + rect.width - margin, yPos, g2D

        GraphicUtils.drawVerticallyCenteredText("Gaps", margin, rect, g2D, true);
        rect.y += rect.height;

        rect.height = EXPANDED_HEIGHT;

        String ref = mgr.getSpeciesName(mgr.refId);
        if (ref == null) {
            ref = mgr.refId;
        }
        GraphicUtils.drawVerticallyCenteredText(ref, margin, rect, g2D, true);
        rect.y += rect.height;

        for (String sp : mgr.getSelectedSpecies()) {

            String name = mgr.getSpeciesName(sp);
            if (name == null) {
                name = sp;
            }

            if (visibleRectangle.intersects(rect)) {

                GraphicUtils.drawVerticallyCenteredText(name, margin, rect, g2D, true);
            }
            rect.y += rect.height;
        }

    }

    public void setColorScale(ContinuousColorScale colorScale) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void render(RenderContext context, Rectangle rect) {

        double locScale = context.getScale();

        if (locScale > 1) {
            Rectangle r = new Rectangle(rect);
            if (visibleNameRect != null) {
                r.y = visibleNameRect.y;
                r.height = visibleNameRect.height;
            }

            Graphics2D g = context.getGraphic2DForColor(Color.black);
            GraphicUtils.drawCenteredText("Zoom in to see alignments.", r, g);
            return;

        }

        double origin = context.getOrigin();
        String chr = context.getChr();

        String mafChr = chrMappings == null ? null : chrMappings.get(chr);
        int start = (int) origin;
        int end = (int) (origin + rect.width * locScale) + 1;
        MAFTile[] tiles = mgr.getTiles(mafChr, start, end);
        if (tiles != null) {
            for (MAFTile tile : tiles) {
                // render tile
                if (tile != null) {
                    renderTile(context, rect, tile);
                }
            }
        }

    }

    private void renderTile(RenderContext context, Rectangle trackRectangle, MAFTile tile) {

        int y = trackRectangle.y;

        MASequence reference = tile.refSeq;
        if (reference == null) {
            return;
        }

        Rectangle rect = new Rectangle(trackRectangle);
        rect.height = GAPS_HEIGHT;
        rect.y = y;
        renderer.renderGaps(tile.getGaps(), context, rect);
        rect.y += rect.height;

        rect.height = EXPANDED_HEIGHT;
        renderer.renderAligment(reference, reference, null, context, rect, this);
        rect.y += rect.height;


        for (String sp : mgr.getSelectedSpecies()) {

            MASequence seq = tile.alignedSequences.get(sp);
            if (seq != null) {
                renderer.renderAligment(seq, reference, tile.getGaps(), context, rect, this);
            }

            rect.y += rect.height;

        }

    }

    public boolean isLogNormalized() {
        return false;
    }

    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {
        return "Multiple alignments";
    }

    public float getRegionScore(String chr, int start, int end, int zoom,
                                RegionScoreType type, String frameName) {
        return 0;
    }


    /**
     * Override to return a specialized popup menu
     *
     * @return
     */
    @Override
    public IGVPopupMenu getPopupMenu(TrackClickEvent te) {
        IGVPopupMenu menu = new IGVPopupMenu();

        JMenuItem configTrack = new JMenuItem("Configure track...");
        configTrack.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent actionEvent) {
                configureTrack();
            }
        });

        menu.add(configTrack);
        return menu;
    }


    @Override
    public boolean handleDataClick(TrackClickEvent te) {
        MouseEvent e = te.getMouseEvent();
        if (e.isPopupTrigger()) {
            configureTrack();
        } else {
            if (IGV.getInstance().isShowDetailsOnClick()) {
                openTooltipWindow(te);
            }
        }
        return true;
    }


    private void configureTrack() {
        MAFConfigurationDialog dialog = new MAFConfigurationDialog(IGV.getMainFrame(), true, mgr);
        dialog.setLocationRelativeTo(IGV.getMainFrame());
        dialog.addWindowListener(new java.awt.event.WindowAdapter() {

            public void windowClosing(java.awt.event.WindowEvent e) {
                System.exit(0);
            }
        });
        dialog.setVisible(true);

        if (dialog.cancelled) {
        } else {
            List<String> selectedSpecies = dialog.getSelectedSpecies();
            mgr.setSelectedSpecies(selectedSpecies);
            PreferenceManager.getInstance().setMafSpecies(selectedSpecies);
            IGV.getInstance().repaint();
        }

    }
}

