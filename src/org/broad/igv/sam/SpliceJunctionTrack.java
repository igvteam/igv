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

import htsjdk.tribble.Feature;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.SpliceJunctionFeature;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.SpliceJunctionRenderer;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.SashimiPlot;
import org.broad.igv.ui.event.AlignmentTrackEvent;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * @author dhmay
 *         Finds splice junctions in real time and renders them as Features
 */
public class SpliceJunctionTrack extends FeatureTrack {

    private static Logger log = Logger.getLogger(SpliceJunctionTrack.class);

    public enum StrandOption {COMBINE, FORWARD, REVERSE, BOTH}

    // Strand option is shared by all tracks
    private static StrandOption strandOption;

    private AlignmentTrack alignmentTrack;
    private IAlignmentDataManager dataManager;
    private boolean removed = false;

    /**
     * The "DataPanel" containing this track.  This field might be null at any given time.  It is updated each repaint.
     */

    private JComponent dataPanel;

    public static void setStrandOption(StrandOption so) {
        strandOption = so;
    }

    public static StrandOption getStrandOption() {
        return strandOption;
    }


    public SpliceJunctionTrack(ResourceLocator locator, String name, IAlignmentDataManager dataManager, AlignmentTrack alignmentTrack, StrandOption ignoreStrand) {
        super(locator, locator.getPath() + "_junctions", name);

        super.setDataRange(new DataRange(0, 0, 60));
        setRendererClass(SpliceJunctionRenderer.class);
        this.dataManager = dataManager;
        this.alignmentTrack = alignmentTrack;
        this.strandOption = ignoreStrand;
        // Register track
    }

    @Override
    protected boolean isShowFeatures(RenderContext context) {
        float maxRange = PreferenceManager.getInstance().getAsFloat(PreferenceManager.SAM_MAX_VISIBLE_RANGE);
        float minVisibleScale = (maxRange * 1000) / 700;
        return context.getScale() < minVisibleScale;
    }

    public boolean isRemoved() {
        return removed;
    }

    public void clear() {
        this.packedFeaturesMap.clear();
    }

    @Override
    public void dispose() {
        super.dispose();
        removed = true;
        dataManager = null;
        alignmentTrack = null;
        setVisible(false);
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


        if (alignmentTrack != null) {
            popupMenu.addSeparator();

            final JMenuItem alignmentItem = new JCheckBoxMenuItem("Show Alignment Track");
            alignmentItem.setSelected(alignmentTrack.isVisible());
            alignmentItem.setEnabled(!alignmentTrack.isRemoved());
            alignmentItem.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    alignmentTrack.onAlignmentTrackEvent(new AlignmentTrackEvent(SpliceJunctionTrack.this, AlignmentTrackEvent.Type.VISIBLE, alignmentItem.isSelected()));
                    if (alignmentItem.isSelected()) {
                        alignmentTrack.onAlignmentTrackEvent(new AlignmentTrackEvent(SpliceJunctionTrack.this, AlignmentTrackEvent.Type.RELOAD));
                    }
                }
            });
            popupMenu.add(alignmentItem);

            final CoverageTrack coverageTrack = alignmentTrack.getCoverageTrack();
            if (coverageTrack != null) {
                final JMenuItem coverageItem = new JCheckBoxMenuItem("Show Coverage Track");
                coverageItem.setSelected(coverageTrack.isVisible());
                coverageItem.setEnabled(!coverageTrack.isRemoved());
                coverageItem.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        UIUtilities.invokeOnEventThread(new Runnable() {

                            public void run() {
                                coverageTrack.setVisible(coverageItem.isSelected());
                                IGV.getInstance().getMainPanel().revalidate();

                            }
                        });
                    }
                });
                popupMenu.add(coverageItem);
            }


            final JMenuItem junctionItem = new JMenuItem("Hide Track");
            junctionItem.setEnabled(!isRemoved());
            junctionItem.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    alignmentTrack.onAlignmentTrackEvent(new AlignmentTrackEvent(SpliceJunctionTrack.this, AlignmentTrackEvent.Type.SPLICE_JUNCTION, false));
                    IGV.getInstance().getMainPanel().revalidate();
                }
            });
            popupMenu.add(junctionItem);

        }


        return popupMenu;
    }

    public boolean isLogNormalized() {
        return false;
    }

    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, String frameName) {
        return 0;
    }

    @Override
    protected String getZoomInMessage(String chr) {
        return "Zoom in to see junctions.";
    }

    @Override
    protected void loadFeatures(String chr, int start, int end, ReferenceFrame referenceFrame) {
        //dataPanel = context.getPanel();
        AlignmentInterval loadedInterval = dataManager.getLoadedInterval(referenceFrame.getCurrentRange());
        if (loadedInterval == null) return;

        SpliceJunctionHelper helper = loadedInterval.getSpliceJunctionHelper();
        List<SpliceJunctionFeature> features = helper.getFilteredJunctions(strandOption);
        if (features == null) {
            features = Collections.emptyList();
        }
        int intervalStart = loadedInterval.getStart();
        int intervalEnd = loadedInterval.getEnd();
        PackedFeatures pf = new PackedFeaturesSpliceJunctions(chr, intervalStart, intervalEnd, features.iterator(), getName());
        packedFeaturesMap.put(referenceFrame.getName(), pf);
      //  if (context.getPanel() != null) context.getPanel().repaint();
    }

    @Override
    public String getExportTrackLine() {
        return "track graphType=junctions";
    }


    @Override
    public boolean handleDataClick(TrackClickEvent te) {
        boolean result = super.handleDataClick(te);
        if (dataPanel != null) dataPanel.repaint();

        return result;
    }

    /**
     * Get all features which overlap the specified locus
     *
     * @return
     */
    @Override
    public List<Feature> getFeatures(String chr, int start, int end) {
        List<Feature> features = new ArrayList<Feature>();
        try {
            Iterator<Feature> iter = source.getFeatures(chr, start, end);
            while (iter.hasNext()) {
                features.add(iter.next());
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return features;
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

                if (dataPanel != null) dataPanel.repaint();
            }
        });

        return autoscaleItem;
    }
    // End of Roche-Tessella modification


}
