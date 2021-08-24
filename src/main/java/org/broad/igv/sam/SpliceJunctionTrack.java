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
import org.broad.igv.Globals;
import org.broad.igv.feature.SpliceJunctionFeature;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.SpliceJunctionRenderer;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.sashimi.SashimiPlot;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.variant.VariantTrack;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * @author dhmay
 * Finds splice junctions in real time and renders them as Features
 */
public class SpliceJunctionTrack extends FeatureTrack {

    private static Logger log = Logger.getLogger(SpliceJunctionTrack.class);

    public enum StrandOption {COMBINE, FORWARD, REVERSE, BOTH}

    // Strand option is shared by all tracks
    private static StrandOption strandOption;

    private AlignmentTrack alignmentTrack;
    private AlignmentDataManager dataManager;
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


    public SpliceJunctionTrack(ResourceLocator locator, String name,
                               AlignmentDataManager dataManager,
                               AlignmentTrack alignmentTrack,
                               StrandOption ignoreStrand) {
        super(locator, locator.getPath() + "_junctions", name);

        super.setDataRange(new DataRange(0, 0, 60));
        setRendererClass(SpliceJunctionRenderer.class);
        if (dataManager != null) {
            dataManager.unsubscribe(this);
        }
        this.dataManager = dataManager;
        this.dataManager.subscribe(this);
        this.alignmentTrack = alignmentTrack;
        this.strandOption = ignoreStrand;
    }

    public SpliceJunctionTrack() {
    }

    @Override
    public String getSample() {
        if (sampleId != null) {
            return sampleId;    // Explicitly set sample ID (e.g. from server load XML)
        }
        return alignmentTrack.getSample();
    }

    protected boolean isShowFeatures(ReferenceFrame frame) {
        float maxRange = PreferencesManager.getPreferences().getAsFloat(Constants.SAM_MAX_VISIBLE_RANGE);
        float minVisibleScale = (maxRange * 1000) / 700;
        return frame.getScale() < minVisibleScale;
    }

    public boolean isRemoved() {
        return removed;
    }

    @Override
    public boolean isVisible() {
        return super.isVisible() && !removed;
    }

    public void clear() {
        this.packedFeaturesMap.clear();
    }

    @Override
    public void dispose() {
        // SpliceJunctionTracks are not disposed of,  merely hidden.
    }

    @Override
    public void setVisible(boolean visible) {
        if (visible != isVisible()) {
            super.setVisible(visible);
            if (visible) {
                dataManager.initLoadOptions();
            }
            if (IGV.hasInstance()) {
                IGV.getInstance().getMainPanel().revalidate();
            }
        }
    }

    /**
     * Override to return a specialized popup menu
     *
     * @return
     */
    @Override
    public IGVPopupMenu getPopupMenu(TrackClickEvent te) {

        IGVPopupMenu menu = new IGVPopupMenu();

        JLabel popupTitle = new JLabel("  " + getName(), JLabel.CENTER);

        Font newFont = menu.getFont().deriveFont(Font.BOLD, 12);
        popupTitle.setFont(newFont);
        if (popupTitle != null) {
            menu.add(popupTitle);
        }
        menu.addSeparator();

        ArrayList<Track> tmp = new ArrayList();
        tmp.add(this);
        TrackMenuUtils.addStandardItems(menu, tmp, te);

        menu.addSeparator();
        menu.add(getChangeAutoScale());


        menu.addSeparator();
        JMenuItem sashimi = new JMenuItem("Sashimi Plot");
        sashimi.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                SashimiPlot.getSashimiPlot(null);
            }
        });
        menu.add(sashimi);


        // Hide/show items
        if (alignmentTrack != null) {

            menu.addSeparator();

            final CoverageTrack coverageTrack = alignmentTrack.getCoverageTrack();
            if (coverageTrack != null) {
                final JMenuItem item = new JCheckBoxMenuItem("Show Coverage Track");
                item.setSelected(coverageTrack.isVisible());
                item.setEnabled(!coverageTrack.isRemoved());
                item.addActionListener(e -> {
                    coverageTrack.setVisible(item.isSelected());
                    IGV.getInstance().repaint(Arrays.asList(coverageTrack));
                });
                menu.add(item);
            }

            final JMenuItem junctionItem = new JCheckBoxMenuItem("Show Splice Junction Track");
            junctionItem.setSelected(true);
            junctionItem.addActionListener(e -> {
                SpliceJunctionTrack.this.setVisible(junctionItem.isSelected());
                IGV.getInstance().repaint(Arrays.asList(SpliceJunctionTrack.this));
            });
            menu.add(junctionItem);
            // Disable if this is the only visible track
            if (!((coverageTrack != null && coverageTrack.isVisible()) ||  alignmentTrack.isVisible())) {
                junctionItem.setEnabled(false);
            }

            final JMenuItem alignmentItem = new JCheckBoxMenuItem("Show Alignment Track");
            alignmentItem.setSelected(alignmentTrack.isVisible());
            alignmentItem.setEnabled(!alignmentTrack.isRemoved());
            alignmentItem.addActionListener(e -> {
                alignmentTrack.setVisible(alignmentItem.isSelected());
                IGV.getInstance().repaint(Arrays.asList(alignmentTrack));
            });
            menu.add(alignmentItem);
        }

        return menu;
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
    protected void renderFeatures(RenderContext context, Rectangle inputRect) {

        // Intercept renderFeatures call and create splice junctions from alignments, if needed.
        ReferenceFrame frame = context.getReferenceFrame();
        if (!packedFeaturesMap.containsKey(frame.getName())) {

            AlignmentInterval loadedInterval = dataManager.getLoadedInterval(frame);
            if (loadedInterval != null) {
                SpliceJunctionHelper helper = loadedInterval.getSpliceJunctionHelper();
                List<SpliceJunctionFeature> features = helper.getFilteredJunctions(strandOption);
                if (features == null) {
                    features = Collections.emptyList();
                }
                int intervalStart = loadedInterval.getStart();
                int intervalEnd = loadedInterval.getEnd();
                PackedFeatures pf = new PackedFeaturesSpliceJunctions(frame.getChrName(), intervalStart, intervalEnd,features.iterator(), getDisplayMode());
                packedFeaturesMap.put(frame.getName(), pf);
            }
        }

        super.renderFeatures(context, inputRect);
    }

    public void load(ReferenceFrame frame) {
        dataManager.load(frame, alignmentTrack.renderOptions, true);

    }

    @Override
    public boolean isReadyToPaint(ReferenceFrame frame) {
        if (frame.getChrName().equals(Globals.CHR_ALL) || frame.getScale() > dataManager.getMinVisibleScale()) {
            return true;   // Nothing to paint
        } else {

            if (!dataManager.isLoaded(frame)) {
                packedFeaturesMap.clear();
                return false;
            } else {
                return true;
            }
        }
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


    @Override
    public void marshalXML(Document document, Element element) {

        super.marshalXML(document, element);

        if (removed) {
            element.setAttribute("removed", String.valueOf(removed));
        }


    }

    @Override
    public void unmarshalXML(Element element, Integer version) {

        super.unmarshalXML(element, version);

        if (element.hasAttribute("removed")) {
            this.removed = Boolean.parseBoolean("removed");
        }

    }

}
