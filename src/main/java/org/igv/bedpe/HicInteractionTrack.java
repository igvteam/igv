package org.igv.bedpe;

import org.igv.hic.HicFile;
import org.igv.renderer.ContinuousColorScale;
import org.igv.track.RenderContext;
import org.igv.track.TrackClickEvent;
import org.igv.ui.panel.FrameManager;
import org.igv.ui.panel.IGVPopupMenu;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.util.ResourceLocator;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import javax.swing.*;
import java.awt.*;
import java.util.List;

/**
 * Subclass of InteractionTrack for HiC format files.
 * Handles HiC-specific behavior like normalization, contact map views, and zoom-based filtering.
 */
public class HicInteractionTrack extends InteractionTrack {

    public HicInteractionTrack() {
        super();
    }

    public HicInteractionTrack(ResourceLocator locator, HicSource source) {
        super(locator, source);
        // HiC-specific defaults
        maxFeatureCount = 5000;
        graphType = GraphType.NESTED_ARC;
        useScore = true;
        setColor(Color.red);
    }

    @Override
    protected List<BedPE> filterFeaturesForZoom(List<BedPE> features, LoadedInterval interval, ReferenceFrame referenceFrame) {
        // In HiC mode we limit interactions to those in view plus a margin of one screen width to either side.
        // If zooming in this means we have to filter the features from the previous zoom level that are outside
        // of this range. Not doing so leads to inconsistent rendering when loading for the current zoom
        // completes and repaints.
        if (interval.zoom() < referenceFrame.getZoom()) {
            int start = (int) referenceFrame.getOrigin();
            int end = (int) referenceFrame.getEnd();
            int w = (end - start);
            int finalStart = start - w;
            int finalEnd = end + w;
            return features.stream()
                    .takeWhile(f -> f.getStart() <= finalEnd)
                    .filter(f -> f.getEnd() >= finalStart)
                    .toList();
        }
        return features;
    }

    @Override
    protected void addFormatSpecificMenuItems(IGVPopupMenu menu, TrackClickEvent te) {
        // Add normalization options for HiC tracks
        List<String> normalizationTypes = featureSource.getNormalizationTypes();
        if (normalizationTypes != null && normalizationTypes.size() > 1) {
            menu.addSeparator();
            menu.add(new JLabel("<html><b>Normalization</b>"));
            ButtonGroup normGroup = new ButtonGroup();
            for (String type : normalizationTypes) {
                String label = normalizationLabels.getOrDefault(type, type);
                JRadioButtonMenuItem normItem = new JRadioButtonMenuItem(label);
                normItem.setSelected(type.equals(normalization));
                normItem.addActionListener(e -> {
                    this.normalization = type;
                    if (contactMapView != null) {
                        contactMapView.setNormalization(type);
                    }
                    this.repaint();
                });
                normGroup.add(normItem);
                menu.add(normItem);
            }
        }

        menu.addSeparator();
        JMenuItem mapItem = new JMenuItem("Open Contact Map View");
        mapItem.setEnabled(contactMapView == null && !FrameManager.isGeneListMode());
        mapItem.addActionListener(e -> {
            ReferenceFrame frame = te.getFrame() != null ? te.getFrame() : FrameManager.getDefaultFrame();
            if (contactMapView == null) {
                ContinuousColorScale colorScale = this.getColorScale();
                HicFile hicFile = ((HicSource) featureSource).getHicFile();
                ContactMapView.showPopup(this, hicFile, normalization, frame, colorScale.getMaxColor());
            }
        });
        menu.add(mapItem);
    }

    @Override
    protected boolean supportsGraphTypeSelection() {
        return false;
    }

    @Override
    protected boolean supportsCircularView() {
        return false;
    }

    @Override
    protected boolean supportsAutoscaleMenu() {
        return false;
    }

    @Override
    protected boolean supportsFeatureWindowMenu() {
        return false;
    }

    @Override
    public void marshalXML(Document document, Element element) {
        super.marshalXML(document, element);

        String nviString = ((HicSource) featureSource).getNVIString();
        if (nviString != null) {
            element.setAttribute("nvi", nviString);
        }
    }

    @Override
    public void unmarshalXML(Element element, Integer version) {
        super.unmarshalXML(element, version);

        if (element.hasAttribute("nvi")) {
            String nviString = element.getAttribute("nvi");
            ((HicSource) featureSource).setNVIString(nviString);
        }
    }
}

