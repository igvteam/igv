/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Fred Hutchinson Cancer Research Center and Broad Institute
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


package org.broad.igv.sam;

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.SpliceJunctionFeature;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.SpliceJunctionRenderer;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.SashimiPlot;
import org.broad.igv.ui.event.AlignmentTrackEvent;
import org.broad.igv.ui.event.AlignmentTrackEventListener;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * @author dhmay
 *         Finds splice junctions in real time and renders them as Features
 */
public class SpliceJunctionFinderTrack extends FeatureTrack {

    private static Logger log = Logger.getLogger(SpliceJunctionFinderTrack.class);

    public enum StrandOption {COMBINE, FORWARD, REVERSE, BOTH}

    // Strand option is shared by all tracks
    private static StrandOption strandOption;

    IAlignmentDataManager dataManager;
    PreferenceManager prefs;

    // The "parent" of the track (a DataPanel).  This release of IGV does not support owner-track relationships
    // directory,  so this field might be null at any given time.  It is updated each repaint.
    JComponent parent;


    public SpliceJunctionFinderTrack(ResourceLocator locator, String name, IAlignmentDataManager dataManager, StrandOption ignoreStrand) {
        super(locator, locator.getPath() + "_junctions", name);

        super.setDataRange(new DataRange(0, 0, 60));
        setRendererClass(SpliceJunctionRenderer.class);
        this.dataManager = dataManager;
        prefs = PreferenceManager.getInstance();
        this.strandOption = ignoreStrand;
        // Register track
    }

    @Override
    protected boolean isShowFeatures(RenderContext context) {
        float maxRange = PreferenceManager.getInstance().getAsFloat(PreferenceManager.SAM_MAX_VISIBLE_RANGE);
        float minVisibleScale = (maxRange * 1000) / 700;
        return context.getScale() < minVisibleScale;
    }


    public static void setStrandOption(StrandOption so) {
        strandOption = so;
    }


    public static StrandOption getStrandOption() {
        return strandOption;
    }

    public void clear() {
        this.packedFeaturesMap.clear();
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

        popupMenu.addSeparator();
        popupMenu.add(getChangeAutoScale());


        popupMenu.addSeparator();
        JMenuItem sashimi = new JMenuItem("Sashimi Plot");
        sashimi.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                SashimiPlot.getSashimiPlot(null);
            }
        });
        popupMenu.add(sashimi);

        return popupMenu;
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
        AlignmentInterval loadedInterval = dataManager.getLoadedInterval(context.getReferenceFrame().getCurrentRange());
        if (loadedInterval == null) return;

        SpliceJunctionHelper helper = loadedInterval.getSpliceJunctionHelper();
        List<SpliceJunctionFeature> features = helper.getFilteredJunctions(strandOption);
        if (features == null) {
            features = Collections.emptyList();
        }
        int intervalStart = loadedInterval.getStart();
        int intervalEnd = loadedInterval.getEnd();
        PackedFeatures pf = new PackedFeaturesSpliceJunctions(chr, intervalStart, intervalEnd, features.iterator(), getName());
        packedFeaturesMap.put(context.getReferenceFrame().getName(), pf);
        if (context.getPanel() != null) context.getPanel().repaint();
    }


    @Override
    public boolean handleDataClick(TrackClickEvent te) {
        boolean result = super.handleDataClick(te);
        if (parent != null) parent.repaint();

        return result;
    }

    // Start of Roche-Tessella modification
    private JMenuItem getChangeAutoScale() {

        final JCheckBoxMenuItem autoscaleItem = new JCheckBoxMenuItem("Autoscale");

        boolean autoScale = getAutoScale();
        autoscaleItem.setSelected(autoScale);


        autoscaleItem.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent evt) {
                boolean autoScale = getAutoScale();
                TrackProperties tp = new TrackProperties();
                if (autoScale) {
                    tp.setAutoScale(false);
                    autoscaleItem.setSelected(false);
                } else {
                    tp.setAutoScale(true);
                    autoscaleItem.setSelected(true);
                }
                tp.setRendererClass(SpliceJunctionRenderer.class);
                setProperties(tp);

                if (parent != null) parent.repaint();
            }
        });

        return autoscaleItem;
    }
    // End of Roche-Tessella modification


}
