package org.igv.sam;

import org.apache.commons.math3.stat.Frequency;
import org.igv.Globals;
import org.igv.feature.IGVFeature;
import org.igv.feature.SpliceJunctionFeature;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.renderer.DataRange;
import org.igv.renderer.GraphicUtils;
import org.igv.renderer.Renderer;
import org.igv.renderer.SpliceJunctionRenderer;
import org.igv.sashimi.SashimiPlot;
import org.igv.track.*;
import org.igv.ui.IGV;
import org.igv.ui.panel.IGVPopupMenu;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.ui.util.MessageUtils;
import org.igv.util.ResourceLocator;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
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
    public TrackType getType() {
        return TrackType.junction;
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

    public boolean isLogNormalized() {
        return false;
    }

    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, String frameName) {
        return 0;
    }

    @Override
    protected void renderFeatures(RenderContext context, Rectangle ignore) {

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
                packedFeaturesMap.put(frame, pf);
            }
        }

        super.renderFeatures(context, ignore);
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
     * Override to return a list of menu items for the popup menu
     *
     * @return
     */
    @Override
    public List<Component> getPopupMenuItems(TrackClickEvent te) {

        List<Component> items = new ArrayList<>();

        JLabel popupTitle = new JLabel("  " + getName(), JLabel.CENTER);
        popupTitle.setFont(UIManager.getFont("Label.font").deriveFont(Font.BOLD, 12));
        items.add(popupTitle);
        items.add(new JPopupMenu.Separator());

        items.add(new JPopupMenu.Separator());
        items.addAll(TrackMenuUtils.getDisplayModeMenuItems(Collections.singletonList(this)));

        items.add(new JPopupMenu.Separator());
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
        items.add(setScaleItem);

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
        items.add(autoscaleItem);

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
        items.add(minJunctionCoverageItem);

        items.add(new JPopupMenu.Separator());
        JMenuItem sashimi = new JMenuItem("Sashimi Plot");
        sashimi.addActionListener(e -> SashimiPlot.openSashimiPlot());
        items.add(sashimi);


        // Hide/show items
        if (alignmentTrack != null) {
            items.add(new JPopupMenu.Separator());
            items.addAll(AlignmentTrackMenuHelper.getShowMenuItems(alignmentTrack));
        }

        return items;
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

    @Override
    public void marshalJSON(org.json.JSONObject json) {
        super.marshalJSON(json);
        if (removed) {
            json.put("removed", removed);
        }

        // Autoscale is set here, but restored in AbstractTrack
        json.put("autoScale", String.valueOf(autoScale));

        json.put("maxdepth", ((SpliceJunctionRenderer) renderer).getMaxDepth());
    }

    @Override
    public void unmarshalJSON(org.json.JSONObject jsonObject) {
        super.unmarshalJSON(jsonObject);
        if (jsonObject.has("removed")) {
            this.removed = jsonObject.getBoolean("removed");
        }
        if (jsonObject.has("maxdepth")) {
            if (renderer != null && renderer instanceof SpliceJunctionRenderer) {
                ((SpliceJunctionRenderer) renderer).setMaxDepth(jsonObject.getInt("maxdepth"));
            }
        }
    }
}
