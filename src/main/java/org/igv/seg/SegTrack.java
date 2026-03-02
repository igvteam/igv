package org.igv.seg;

import org.igv.feature.LocusScore;
import org.igv.feature.genome.Genome;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.renderer.*;
import org.igv.renderer.Renderer;
import org.igv.sample.SampleGroup;
import org.igv.sample.SampleMenuUtils;
import org.igv.session.RendererFactory;
import org.igv.track.*;
import org.igv.ui.FontManager;
import org.igv.ui.panel.IGVPopupMenu;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.util.ResourceLocator;
import org.json.JSONObject;
import org.w3c.dom.Element;

import java.awt.*;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.*;
import java.util.List;

public class SegTrack extends AbstractTrack {

    private static final Logger log = LogManager.getLogger(SegTrack.class);

    SegmentedDataSet dataset;
    int sampleHeight = 15;
    Renderer<LocusScore> renderer = new HeatmapRenderer();

    public SegTrack(ResourceLocator locator, String id, String name, SegmentedDataSet dataset, Genome genome) {

        super(locator, id, name);
        this.dataset = dataset;
        setDataRange(new DataRange(0, 1, 2));
        setDataType(dataset.getType());
        initScale(dataset);

        this.initSamples(dataset.getSampleNames());
    }

    public TrackType getType() {
        return TrackType.seg;
    }

    void initScale(SegmentedDataSet dataset) {

        var min = (float) dataset.getDataMin();
        var max = (float) dataset.getDataMax();
        var baseline = 0;

        // If the range is all positive numbers set the min to zero
        if (min > 0) {
            min = 0;
        }

        setDataRange(new DataRange(min, baseline, max));
    }


    @Override
    public int getContentHeight() {
        var nSamples = sampleCount();
        return nSamples * sampleHeight + (getSampleGroups().size() - 1) * groupGap;
    }

    @Override
    public int getSampleHeight() {
        return sampleHeight;
    }

    @Override
    public void render(RenderContext context) {

        Rectangle clipBounds = context.getClipBounds();
        Rectangle trackRect = context.getTrackRectangle();
        var y = 0;
        for (var group : getSampleGroups()) {
            for (String sample : group.samples()) {
                if (y > clipBounds.y + clipBounds.height) {
                    break;
                }
                if (y + sampleHeight > clipBounds.y) {
                    var r = new Rectangle(trackRect.x, y, trackRect.width, sampleHeight);
                    var scores = dataset.getSegments(sample, context.getChr());
                    this.renderer.render(scores, context, r, this);
                }
                y += sampleHeight;
            }
            y += groupGap; // Gap between groups
        }
    }

    @Override
    public void renderName(Graphics2D g2D, Rectangle visibleRect, Rectangle clipRect) {

        // Calculate fontsize
        var gap = Math.min(4, sampleHeight / 3);
        var fs = Math.min(fontSize, sampleHeight - gap);
        var font = FontManager.getFont(fs);
        g2D.setFont(font);

        var y = 0;
        for (var group : getSampleGroups()) {
            for (String s : group.samples()) {
                if (y > clipRect.y + clipRect.height) {
                    break;
                }
                if (y + sampleHeight > clipRect.y) {
                    var r = new Rectangle(0, y, visibleRect.width, sampleHeight);
                    GraphicUtils.drawWrappedText(s, r, g2D, false);
                }
                y += sampleHeight;
            }
            y += groupGap; // Gap between groups
        }
    }

    @Override
    public void setRendererClass(Class rc) {
        try {
            if (Renderer.class.isAssignableFrom(rc)) {
                renderer = (Renderer<LocusScore>) rc.getDeclaredConstructor().newInstance();
            }
        } catch (Exception ex) {
            log.error("Error instantiating renderer ", ex);
        }
    }

    public IGVPopupMenu getPopupMenu(final TrackClickEvent te) {
        var menu = TrackMenuUtils.getPopupMenu(Collections.singletonList(this), "Menu", te);
        menu.addSeparator();
        TrackMenuUtils.addDataRendererItems(menu, Arrays.asList("Heatmap", "Bar Chart", "Points", "Line Plot"), Collections.singletonList(this));


        List<String> keys = AttributeManager.getInstance().getAttributeNames();

        if(keys.size() > 0) {
            menu.addSeparator();
            menu.add(SampleMenuUtils.getSortByAttributeItem(this));
            menu.add(SampleMenuUtils.getGroupByAttributeItem(this));
            menu.add(SampleMenuUtils.getFilterByAttributeItem(this));
        }

        return menu;
    }

    /**
     * Return a value string for the tooltip window at the given location, or null to signal there is no value
     * at that location
     *
     * @param chr      The chromosome
     * @param position The genomic position
     * @param mouseX   The mouse X coordinate
     * @param mouseY   Mouse position relative to the enclosing panel
     * @param frame    The reference frame
     * @return A value string for the tooltip, or null.
     */
    @Override
    public String getValueStringAt(String chr, double position, int mouseX, int mouseY, ReferenceFrame frame) {

        var trackY = mouseY - this.getY();
        var y = 0;
        for (var group : getSampleGroups()) {
            var groupPixelHeight = group.samples().size() * sampleHeight;
            if (trackY >= y && trackY < y + groupPixelHeight) {
                var sampleIndex = (trackY - y) / sampleHeight;
                var sampleName = group.samples().get(sampleIndex);
                var buf = new StringBuilder();
                var score = dataset.getSegmentAt(sampleName, chr, position, frame);
                // If there is no value here, return null to signal no popup
                if (score == null) {
                    return null;
                }
                buf.append(sampleName).append("<br>");
                if ((getDataRange() != null) && (getRenderer() instanceof XYPlotRenderer)) {
                    buf.append("Data scale: ").append(getDataRange().getMinimum()).append(" - ").append(getDataRange().getMaximum()).append("<br>");
                }

                buf.append(score.getValueString(position, mouseX, getWindowFunction()));
                return buf.toString();
            } else if (trackY < y + groupPixelHeight + groupGap) {
                return null; // In a gap
            }
            y += groupPixelHeight + groupGap;
        }
        return null;
    }


    /**
     * Get the score over the provided region for the given type. Different types
     * are processed differently. Results are cached according to the provided frameName,
     * if provided. If not, a string is created based on the inputs.
     *
     * @param chr   The chromosome
     * @param start The start of the genomic range
     * @param end   The end of the genomic range
     * @param type  The score type
     */
    @Override
    public void sortSamplesByValue(String chr, int start, int end, RegionScoreType type) {

        if (end <= start) {
            return;
        }

        for (var group : getSampleGroups()) {
            // Compute a value for each sample
            var sampleScores = new HashMap<String, Float>();
            for (String s : group.samples()) {
                var scores = dataset.getSegments(s, chr);
                var regionScore = 0f;
                var intervalSum = 0;
                var hasNan = false;
                for (var score : scores) {
                    if ((score.getEnd() >= start) && (score.getStart() <= end)) {
                        var interval = Math.min(end, score.getEnd()) - Math.max(start, score.getStart());
                        var value = score.getScore();
                        //For sorting it makes sense to skip NaNs. Not sure about other contexts
                        if (Float.isNaN(value)) {
                            hasNan = true;
                            continue;
                        }
                        regionScore += value * interval;
                        intervalSum += interval;
                    }
                }
                if (intervalSum <= 0) {
                    if (hasNan) {
                        //If the only existing scores are NaN, the overall score should be NaN
                        sampleScores.put(s, Float.NaN);
                    } else {
                        // No scores in interval
                        sampleScores.put(s, -Float.MAX_VALUE);
                    }
                } else {
                    regionScore /= intervalSum;
                    sampleScores.put(s, (type == RegionScoreType.DELETION) ? -regionScore : regionScore);
                }
            }

            group.samples().sort((o1, o2) -> {
                var v1 = sampleScores.get(o1);
                var v2 = sampleScores.get(o2);
                // Put NaNs at the end
                if (v1.isNaN()) {
                    return (v2.isNaN()) ? 0 : 1;
                } else if (v2.isNaN()) {
                    return -1;
                }
                return -Float.compare(v1, v2);
            });
        }
    }

    /**
     * The track itself is not filterable, but samples are.
     *
     * @return false
     */
    @Override
    public boolean isFilterable() {
        return false;
    }

    public static boolean isSegFile(String pathOrURL) {
        String path;
        try {
            var url = new URL(pathOrURL);
            path = url.getPath();
        } catch (MalformedURLException e) {
            // Not a valid URL, assume it's a file path
            path = pathOrURL;
        }
        var lowerPath = path.toLowerCase();
        return lowerPath.endsWith(".seg") || lowerPath.endsWith(".seg.gz");
    }

    @Override
    public void unmarshalXML(Element element, Integer version) {

        super.unmarshalXML(element, version);
        if (element.hasAttribute("renderer")) {

            Class rendererClass = RendererFactory.getRendererClass(element.getAttribute("renderer"));
            if (rendererClass != null) {
                try {
                    renderer = (DataRenderer) rendererClass.newInstance();

                } catch (Exception e) {
                    log.error("Error instantiating renderer: " + rendererClass.getName(), e);
                }
            }
        }
    }

    @Override
    public void marshalJSON(JSONObject json) {
        super.marshalJSON(json);
        if (renderer != null) {
            json.put("renderer", renderer.getClass().getName());
        }
    }

    @Override
    public void unmarshalJSON(JSONObject jsonObject) {
        super.unmarshalJSON(jsonObject);
        if (jsonObject.has("renderer")) {
            Class rendererClass = RendererFactory.getRendererClass(jsonObject.getString("renderer"));
            if (rendererClass != null) {
                try {
                    renderer = (DataRenderer) rendererClass.newInstance();
                } catch (Exception e) {
                    log.error("Error instantiating renderer: " + rendererClass.getName(), e);
                }
            }
        }
    }
}
