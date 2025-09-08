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

import org.apache.commons.math3.stat.Frequency;
import org.broad.igv.Globals;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.SpliceJunctionFeature;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.renderer.SpliceJunctionRenderer;
import org.broad.igv.sashimi.SashimiPlot;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ResourceLocator;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * @author dhmay
 * Finds splice junctions in real time and renders them as Features
 */
public class SpliceJunctionTrack extends FeatureTrack implements ScalableTrack {

    private static Logger log = LogManager.getLogger(SpliceJunctionTrack.class);

    public enum StrandOption {COMBINE, FORWARD, REVERSE, BOTH}

    // Strand option is shared by all tracks
    private static StrandOption strandOption;

    public static void setStrandOption(StrandOption so) {
        strandOption = so;
    }

    public static StrandOption getStrandOption() {
        return strandOption;
    }

    private AlignmentTrack alignmentTrack;
    private AlignmentDataManager dataManager;
    private boolean removed = false;

    public SpliceJunctionTrack(ResourceLocator locator, String name,
                               AlignmentDataManager dataManager,
                               AlignmentTrack alignmentTrack,
                               StrandOption ignoreStrand) {
        super(locator, locator.getPath() + "_junctions", name);

        super.setDataRange(new DataRange(0, 0, 60));
        this.renderer = new SpliceJunctionRenderer();
        if (dataManager != null) {
            dataManager.unsubscribe(this);
        }
        this.dataManager = dataManager;
        this.dataManager.subscribe(this);
        this.alignmentTrack = alignmentTrack;
        SpliceJunctionTrack.strandOption = ignoreStrand;
    }

    public SpliceJunctionTrack() {
        this.renderer = new SpliceJunctionRenderer();
    }

    @Override
    public void setRenderer(Renderer renderer) {
        this.renderer = renderer;
    }

    @Override
    public String getSample() {
        if (getSampleId() != null) {
            return getSampleId();    // Explicitly set sample ID (e.g. from server load XML)
        }
        return alignmentTrack.getSample();
    }

    protected boolean isShowFeatures(ReferenceFrame frame) {
        return frame.getCurrentRange().getLength() <= dataManager.getVisibilityWindow();
    }

    @Override
    public void render(RenderContext context, Rectangle rect) {
        if (!isShowFeatures(context.getReferenceFrame())) {
            Rectangle visibleRect = context.getVisibleRect().intersection(rect);
            Graphics2D g = context.getGraphic2DForColor(Color.gray);
            String message = context.getReferenceFrame().getChrName().equals(Globals.CHR_ALL) ?
                    "Select a chromosome and zoom in to see features." :
                    "Zoom in to see features.";
            GraphicUtils.drawCenteredText(message, visibleRect, g);
        } else {
            super.render(context, rect);
        }
    }


    /**
     * The ScalableTrack interface.  Return a Range object reflecting min/max of the data (splice junction counts
     * in this case.  Used for autoscaling.
     *
     * @param referenceFrame
     * @return
     */
    @Override
    public Range getInViewRange(ReferenceFrame referenceFrame) {

        PackedFeatures packedFeatures = packedFeaturesMap.get(referenceFrame.getName());
        if (packedFeatures == null) {
            return new Range(0, ((SpliceJunctionRenderer) renderer).getMaxDepth());
        } else {

            Frequency f = new Frequency();
            List<Integer> scores = new ArrayList<Integer>();

            List<IGVFeature> featureList = packedFeatures.getFeatures();
            for (IGVFeature feature : featureList) {
                SpliceJunctionFeature junctionFeature = (SpliceJunctionFeature) feature;
                if (feature.getEnd() >= referenceFrame.getOrigin() && feature.getStart() <= referenceFrame.getEnd()) {
                    f.addValue(junctionFeature.getScore());
                    scores.add((int) junctionFeature.getScore());
                }
            }

            Collections.sort(scores);
            Collections.reverse(scores);
            int max = 0;
            for (int s : scores) {
                if (f.getCumPct(s) < 0.99) {
                    max = s;
                    break;
                }
            }
            return new Range(0, max);
        }
    }

    @Override
    public void setDataRange(DataRange axisDefinition) {
        ((SpliceJunctionRenderer) renderer).setMaxDepth((int) axisDefinition.getMaximum());
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
    public void setVisible(boolean visible) {
        if (visible != isVisible()) {
            super.setVisible(visible);
            if (IGV.hasInstance()) {
                IGV.getInstance().getMainPanel().revalidate();
            }
        }
    }

    public boolean isLogNormalized() {
        return false;
    }

    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, String frameName) {
        return 0;
    }

    @Override
    protected void renderFeatures(RenderContext context, Rectangle inputRect) {

        // Intercept renderFeatures call and create splice junctions from alignments, if needed.
        ReferenceFrame frame = context.getReferenceFrame();
        if (!packedFeaturesMap.containsKey(frame.getName())) {

            AlignmentInterval loadedInterval = dataManager.getLoadedInterval(frame);
            if (loadedInterval != null) {
                SpliceJunctionHelper helper = loadedInterval.getSpliceJunctionHelper();
                int minJunctionCoverage = alignmentTrack.getRenderOptions().getMinJunctionCoverage();
                List<SpliceJunctionFeature> features = helper.getFilteredJunctions(strandOption, minJunctionCoverage);

                if (features == null) {
                    features = Collections.emptyList();
                }
                int intervalStart = loadedInterval.getStart();
                int intervalEnd = loadedInterval.getEnd();
                PackedFeatures pf = new PackedFeaturesSpliceJunctions(frame.getChrName(), intervalStart, intervalEnd, features.iterator(), getDisplayMode());
                packedFeaturesMap.put(frame.getName(), pf);
            }
        }

        super.renderFeatures(context, inputRect);
    }

    public void load(ReferenceFrame frame) {
        dataManager.load(frame, alignmentTrack.getRenderOptions(), alignmentTrack.getDisplayMode(), true);

    }

    @Override
    public boolean isReadyToPaint(ReferenceFrame frame) {
        double extent = frame.getEnd() - frame.getOrigin();
        if (frame.getChrName().equals(Globals.CHR_ALL) || extent > dataManager.getVisibilityWindow()) {
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


    public int getMinJunctionCoverage() {
        return alignmentTrack.getRenderOptions().getMinJunctionCoverage();
    }

    public void setMinJunctionCoverage(int minJunctionCoverage) {
        alignmentTrack.getRenderOptions().setMinJunctionCoverage(minJunctionCoverage);
        this.packedFeaturesMap.clear();
        this.repaint();
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

        TrackMenuUtils.addSharedItems(menu, Collections.singletonList(this));

        menu.addSeparator();
        TrackMenuUtils.addDisplayModeItems(Collections.singletonList(this), menu);

        menu.addSeparator();
        final JMenuItem setScaleItem = new JMenuItem("Set Maximum Depth");
        final JCheckBoxMenuItem autoscaleItem = new JCheckBoxMenuItem("Autoscale");

        setScaleItem.addActionListener(evt -> {
            Integer newDepth = TrackMenuUtils.getIntegerInput("Maximum Depth", ((SpliceJunctionRenderer) renderer).getMaxDepth());
            if (newDepth != null && newDepth > 0) {
                ((SpliceJunctionRenderer) renderer).setMaxDepth(newDepth);
                setAutoScale(false);
                autoscaleItem.setSelected(false);
                this.repaint();
            }
        });
        menu.add(setScaleItem);

        autoscaleItem.setSelected(getAutoScale());
        autoscaleItem.addActionListener(evt -> {
            if (getAutoScale()) {
                setAutoScale(false);
                autoscaleItem.setSelected(false);
            } else {
                setAutoScale(true);
                autoscaleItem.setSelected(true);
            }
            this.repaint();
        });
        menu.add(autoscaleItem);

        JMenuItem minJunctionCoverageItem = new JMenuItem("Set Junction Coverage Min");
        minJunctionCoverageItem.setToolTipText("Junctions below this threshold will be removed from view");
        minJunctionCoverageItem.addActionListener(e1 -> {
            int minCov = alignmentTrack.getRenderOptions().getMinJunctionCoverage();
            String input = JOptionPane.showInputDialog("Set Minimum Junction Coverage", minCov);
            if (input == null || input.length() == 0) return;
            try {
                int newMinJunctionCoverage = Integer.parseInt(input);
                setMinJunctionCoverage(newMinJunctionCoverage);
            } catch (NumberFormatException ex) {
                MessageUtils.showMessage("" + input + " is not an integer");
            }
        });
        menu.add(minJunctionCoverageItem);

        menu.addSeparator();
        JMenuItem sashimi = new JMenuItem("Sashimi Plot");
        sashimi.addActionListener(e -> SashimiPlot.openSashimiPlot());
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
            if (!((coverageTrack != null && coverageTrack.isVisible()) || alignmentTrack.isVisible())) {
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


    @Override
    public void marshalXML(Document document, Element element) {

        super.marshalXML(document, element);

        if (removed) {
            element.setAttribute("removed", String.valueOf(removed));
        }

        // Autoscale is set here, but restored in AbstractTrack
        element.setAttribute("autoScale", String.valueOf(autoScale));

        element.setAttribute("maxdepth", String.valueOf(((SpliceJunctionRenderer) renderer).getMaxDepth()));


    }

    @Override
    public void unmarshalXML(Element element, Integer version) {

        super.unmarshalXML(element, version);

        if (element.hasAttribute("removed")) {
            this.removed = Boolean.parseBoolean(element.getAttribute("removed"));
        }

        if (element.hasAttribute("maxdepth")) {
            if (renderer != null && renderer instanceof SpliceJunctionRenderer) {
                ((SpliceJunctionRenderer) renderer).setMaxDepth((int) Float.parseFloat(element.getAttribute("maxdepth")));
            }
        }

    }

}
