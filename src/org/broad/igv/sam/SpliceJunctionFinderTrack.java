/**
 * Copyright (c) 2010-2011 by Fred Hutchinson Cancer Research Center.  All Rights Reserved.

 * This software is licensed under the terms of the GNU Lesser General
 * Public License (LGPL), Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.

 * THE SOFTWARE IS PROVIDED "AS IS." FRED HUTCHINSON CANCER RESEARCH CENTER MAKES NO
 * REPRESENTATIONS OR WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED,
 * INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 * PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS,
 * WHETHER OR NOT DISCOVERABLE.  IN NO EVENT SHALL FRED HUTCHINSON CANCER RESEARCH
 * CENTER OR ITS TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR
 * ANY DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR
 * CONSEQUENTIAL DAMAGES, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS,
 * REGARDLESS OF  WHETHER FRED HUTCHINSON CANCER RESEARCH CENTER SHALL BE ADVISED,
 * SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE
 * FOREGOING.
 */

/*
 *
 */
package org.broad.igv.sam;

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.SpliceJunctionFeature;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.SpliceJunctionRenderer;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.event.AlignmentTrackEvent;
import org.broad.igv.ui.event.AlignmentTrackEventListener;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

/**
 * @author dhmay
 *         Finds splice junctions in real time and renders them as Features
 */
public class SpliceJunctionFinderTrack extends FeatureTrack implements AlignmentTrackEventListener {

    private static Logger log = Logger.getLogger(SpliceJunctionFinderTrack.class);

    AlignmentDataManager dataManager;
    PreferenceManager prefs;
    RenderContext context;
    Genome genome;

    // The "parent" of the track (a DataPanel).  This release of IGV does not support owner-track releationships
    // directory,  so this field might be null at any given time.  It is updated each repaint.
    JComponent parent;


    public SpliceJunctionFinderTrack(ResourceLocator locator, String name, AlignmentDataManager dataManager, Genome genome) {
        super(locator, locator.getPath() + "_junctions", name);

        super.setDataRange(new DataRange(0, 0, 60));
        this.genome = genome;
        setRendererClass(SpliceJunctionRenderer.class);
        this.dataManager = dataManager;
        prefs = PreferenceManager.getInstance();
        // Register track
        IGV.getInstance().addAlignmentTrackEventListener(this);
    }


    @Override
    protected boolean isShowFeatures(RenderContext context) {
        float maxRange = PreferenceManager.getInstance().getAsFloat(PreferenceManager.SAM_MAX_VISIBLE_RANGE);
        float minVisibleScale = (maxRange * 1000) / 700;
        return context.getScale() < minVisibleScale;

    }


    /**
     * Override to return a specialized popup menu
     *
     * @return
     */
    @Override
    public IGVPopupMenu getPopupMenu(TrackClickEvent te) {

        IGVPopupMenu popupMenu = new IGVPopupMenu();

        JLabel popupTitle = new JLabel("  " + getName(), JLabel.CENTER);

        Font newFont = popupMenu.getFont().deriveFont(Font.BOLD, 12);
        popupTitle.setFont(newFont);
        if (popupTitle != null) {
            popupMenu.add(popupTitle);
        }

        popupMenu.addSeparator();

        ArrayList<Track> tmp = new ArrayList();
        tmp.add(this);
        TrackMenuUtils.addStandardItems(popupMenu, tmp, te);
        return popupMenu;
    }


    @Override
    public void setDataRange(DataRange axisDefinition) {
        // Explicitly setting a data range turns off auto-scale
        super.setDataRange(axisDefinition);
    }


    public boolean isLogNormalized() {
        return false;
    }

    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, String frameName) {
        return 0;
    }

    @Override
    protected void loadFeatures(String chr, int start, int end, RenderContext context) {
        parent = context.getPanel();
        final Collection<AlignmentInterval> loadedIntervals = dataManager.getLoadedIntervals(context.getReferenceFrame());
        if(loadedIntervals == null) return;

        for (AlignmentInterval loadedInterval : loadedIntervals) {
            if (loadedInterval != null) {
                List<SpliceJunctionFeature> features = loadedInterval.getSpliceJunctions();
                if (features == null) {
                    features = Collections.emptyList();
                }
                int intervalStart = loadedInterval.getStart();
                int intervalEnd = loadedInterval.getEnd();
                PackedFeatures pf = new PackedFeaturesSpliceJunctions(chr, intervalStart, intervalEnd, features.iterator(), getName());
                packedFeaturesMap.put(context.getReferenceFrame().getName(), pf);
                if (context.getPanel() != null) context.getPanel().repaint();
            }
        }
    }


    @Override
    public boolean handleDataClick(TrackClickEvent te) {
        boolean result = super.handleDataClick(te);
        if (parent != null) parent.repaint();

        return result;
    }

    public void onAlignmentTrackEvent(AlignmentTrackEvent e) {
        AlignmentTrackEvent.Type type = e.getType();
        switch (type) {
            case SPLICE_JUNCTION:
                packedFeaturesMap.clear();
        }

    }
}
