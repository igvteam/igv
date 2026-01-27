package org.igv.seg;

import org.igv.feature.LocusScore;
import org.igv.feature.genome.Genome;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.renderer.*;
import org.igv.session.RendererFactory;
import org.igv.track.*;
import org.igv.ui.FontManager;
import org.igv.ui.IGV;
import org.igv.ui.panel.AttributeHeaderPanel;
import org.igv.ui.panel.IGVPopupMenu;
import org.igv.ui.panel.MouseableRegion;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.util.ResourceLocator;
import org.igv.util.TrackFilter;
import org.w3c.dom.Document;
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
    int groupGap = 15;
    List<SampleGroup> sampleGroups;
    Renderer<LocusScore> renderer = new HeatmapRenderer();


    public SegTrack(ResourceLocator locator, String id, String name, SegmentedDataSet dataset, Genome genome) {
        super(locator, id, name);
        this.dataset = dataset;
        var sampleNames = new ArrayList<>(dataset.getSampleNames());
        this.sampleGroups = new ArrayList<>();
        this.sampleGroups.add(new SampleGroup("", sampleNames));
        setDataRange(new DataRange(0, 1, 2));
        setTrackType(dataset.getType());
        initScale(dataset);

        // Groups by global attribute if it exists
        groupSamplesByAttribute(null);
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
        var nSamples = 0;
        for (var group : sampleGroups) {
            nSamples += group.samples().size();
        }
        return nSamples * sampleHeight + (sampleGroups.size() - 1) * groupGap;
    }

    @Override
    public  void render(RenderContext context, Rectangle visibleRect) {
System.out.println("SegTrack.render: " + getName() + " " + visibleRect);
        var y = 0;
        for (var group : sampleGroups) {
            for (String sample : group.samples()) {
                if(y > visibleRect.y + visibleRect.height) {
                    break;
                }
                if(y + sampleHeight > visibleRect.y) {
                    var r = new Rectangle(visibleRect.x, y, visibleRect.width, sampleHeight);
                    var scores = dataset.getSegments(sample, context.getChr());
                    this.renderer.render(scores, context, r, this);
                }
                y += sampleHeight;
            }
            y += groupGap; // Gap between groups
        }
    }

    @Override
    public void renderAttributes(Graphics2D graphics, Rectangle trackRectangle, Rectangle visibleRect,
                                 List<String> attributeNames, List<MouseableRegion> mouseRegions) {

        final var attributeManager = AttributeManager.getInstance();

        var y = trackRectangle.y;
        for (var group : sampleGroups) {
            for (String s : group.samples()) {
                if (y + sampleHeight > trackRectangle.y && y < trackRectangle.y + trackRectangle.height) {
                    var x = trackRectangle.x;
                    for (var name : attributeNames) {
                        final var key = name.toUpperCase();
                        var attributeValue = attributeManager.getAttribute(s, key);
                        if (attributeValue != null) {
                            var rect = new Rectangle(x, y, AttributeHeaderPanel.ATTRIBUTE_COLUMN_WIDTH, sampleHeight);
                            graphics.setColor(AttributeManager.getInstance().getColor(key, attributeValue));
                            graphics.fill(rect);
                            mouseRegions.add(new MouseableRegion(rect, key, attributeValue));
                        }
                        x += AttributeHeaderPanel.ATTRIBUTE_COLUMN_WIDTH + AttributeHeaderPanel.COLUMN_BORDER_WIDTH;
                    }
                }
                y += sampleHeight;
            }
            y += groupGap; // Gap
        }
    }

    @Override
    public void renderName(Graphics2D g2D, Rectangle trackRectangle) {

        // Calculate fontsize
        var gap = Math.min(4, sampleHeight / 3);
        var fs = Math.min(fontSize, sampleHeight - gap);
        var font = FontManager.getFont(fs);
        g2D.setFont(font);

        var y = 0;
        for (var group : sampleGroups) {
            for (String s : group.samples()) {
                if(y > trackRectangle.y + trackRectangle.height) {
                    break;
                }
                if(y + sampleHeight > trackRectangle.y) {
                    var r = new Rectangle(0, y, trackRectangle.width, sampleHeight);
                    GraphicUtils.drawWrappedText(s, r, g2D, false);
                }
                y += sampleHeight;
            }
            y += groupGap; // Gap between groups
        }
    }


    /**
     * Sort samples.  Sort both the master list and groups, if any.
     *
     * @param comparator the comparator to sort by
     */
    public void sortSamplesByAttribute(Comparator<String> comparator) {
        for (var group : sampleGroups) {
            group.samples().sort(comparator);
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
        for (var group : sampleGroups) {
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

    @Override
    public int sampleCount() {
        var count = 0;
        for (var group : sampleGroups) {
            count += group.samples().size();
        }
        return count;
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

        for (var group : sampleGroups) {
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

    @Override
    public void filterSamples(TrackFilter trackFilter) {
        if (trackFilter == null || trackFilter.isShowAll()) {
            var sampleNames = new ArrayList<>(dataset.getSampleNames());
            this.sampleGroups.clear();
            this.sampleGroups.add(new SampleGroup("All Samples", sampleNames));
        } else {
            var filteredSamples = trackFilter.evaluateSamples(dataset.getSampleNames());
            this.sampleGroups.clear();
            this.sampleGroups.add(new SampleGroup("Filtered Samples", filteredSamples));
        }
    }

    @Override
    public void groupSamplesByAttribute(String attributeKey) {

        if(attributeKey == null && IGV.hasInstance()) {
            attributeKey = IGV.getInstance().getGroupByAttribute();
        }

        this.sampleGroups.clear();

        if (attributeKey == null) {
            var sampleNames = new ArrayList<>(dataset.getSampleNames());
            this.sampleGroups.add(new SampleGroup("", sampleNames));
            repaint();
        } else {
            var attributeManager = AttributeManager.getInstance();
            var groupMap = new LinkedHashMap<String, List<String>>();
            for (var sample : dataset.getSampleNames()) {
                var attributeValue = attributeManager.getAttribute(sample, attributeKey.toUpperCase());
                if (attributeValue == null) {
                    attributeValue = "";
                }
                var samples = groupMap.computeIfAbsent(attributeValue, k -> new ArrayList<>());
                samples.add(sample);
            }

            // Create groups
            for (var key : groupMap.keySet()) {
                sampleGroups.add(new SampleGroup(key, groupMap.get(key)));
            }
        }

        repaint();
    }

    public void marshalXML(Document document, Element element) {
        super.marshalXML(document, element);
        if (renderer != null) {
            var type = RendererFactory.getRenderType(renderer);
            if (type != null) {
                element.setAttribute("renderer", type.name());
            }
        }
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

}
