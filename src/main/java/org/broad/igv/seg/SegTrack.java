package org.broad.igv.seg;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.renderer.*;
import org.broad.igv.track.*;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.panel.AttributeHeaderPanel;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.MouseableRegion;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TrackFilter;

import java.awt.*;
import java.util.*;
import java.util.List;

public class SegTrack extends AbstractTrack {

    private static Logger log = LogManager.getLogger(SegTrack.class);

    SegmentedDataSet dataset;
    Color color = Color.RED;
    Color altColor = Color.BLUE;
    int sampleHeight = 15;
    List<String> sampleNames;
    Renderer renderer = new HeatmapRenderer();


    public SegTrack(ResourceLocator locator, String id, String name, SegmentedDataSet dataset, Genome genome) {
        super(locator, id, name);
        this.dataset = dataset;
        this.sampleNames = new ArrayList(dataset.getSampleNames());
        setDataRange(new DataRange(0, 1, 2));
        setTrackType(dataset.getType());
        initScale(dataset);
    }

    void initScale(SegmentedDataSet dataset) {

        float min = (float) dataset.getDataMin();
        float max = (float) dataset.getDataMax();
        float baseline = 0;

        // If the range is all positive numbers set the min to zero
        if (min > 0) {
            min = 0;
        }

        setDataRange(new DataRange(min, baseline, max));
    }

    @Override
    public int getHeight() {
        return sampleNames.size() * sampleHeight;
    }

    @Override
    public boolean isReadyToPaint(ReferenceFrame frame) {
        return true;
    }

    @Override
    public void render(RenderContext context, Rectangle rect) {

        int y = rect.y;
        for (String s : sampleNames) {
            Rectangle r = new Rectangle(rect.x, y, rect.width, sampleHeight);
            List<LocusScore> scores = dataset.getSegments(s, context.getChr());
            this.renderer.render(scores, context, r, this);
            y += sampleHeight;
        }
    }

    @Override
    public void renderAttributes(Graphics2D graphics, Rectangle trackRectangle, Rectangle visibleRect,
                                 List<String> attributeNames, List<MouseableRegion> mouseRegions) {

        final AttributeManager attributeManager = AttributeManager.getInstance();

        int y = trackRectangle.y;
        for (String s : sampleNames) {
            if (y + sampleHeight > trackRectangle.y && y < trackRectangle.y + trackRectangle.height) {
                int x = trackRectangle.x;
                for (String name : attributeNames) {
                    final String key = name.toUpperCase();
                    String attributeValue = attributeManager.getAttribute(s, key);
                    if (attributeValue != null) {
                        Rectangle rect = new Rectangle(x, y, AttributeHeaderPanel.ATTRIBUTE_COLUMN_WIDTH, sampleHeight);
                        graphics.setColor(AttributeManager.getInstance().getColor(key, attributeValue));
                        graphics.fill(rect);
                        mouseRegions.add(new MouseableRegion(rect, key, attributeValue));
                    }
                    x += AttributeHeaderPanel.ATTRIBUTE_COLUMN_WIDTH + AttributeHeaderPanel.COLUMN_BORDER_WIDTH;
                }
            }
            y += sampleHeight;
        }
    }

    @Override
    public void renderName(Graphics2D g2D, Rectangle trackRectangle, Rectangle visibleRectangle) {

        // Calculate fontsize
        int gap = Math.min(4, sampleHeight / 3);
        int fs = Math.min(fontSize, sampleHeight - gap);
        Font font = FontManager.getFont(fs);
        g2D.setFont(font);

        int y = trackRectangle.y;
        for (String s : sampleNames) {
            Rectangle r = new Rectangle(trackRectangle.x, y, trackRectangle.width, sampleHeight);
            GraphicUtils.drawWrappedText(s, r, g2D, false);
            y += sampleHeight;
        }
    }

    @Override
    public boolean hasSamples() {
        return super.hasSamples();
    }

    /**
     * Sort samples.  Sort both the master list and groups, if any.
     *
     * @param comparator the comparator to sort by
     */
    public void sortSamplesByAttribute(Comparator<String> comparator) {
        if (sampleNames != null) {
            Collections.sort(sampleNames, comparator);
        }
    }

    @Override
    public void setRendererClass(Class rc) {
        try {
            renderer = (DataRenderer) rc.getDeclaredConstructor().newInstance();
        } catch (Exception ex) {
            log.error("Error instantiating renderer ", ex);
        }
    }

    public IGVPopupMenu getPopupMenu(final TrackClickEvent te) {
        IGVPopupMenu menu = TrackMenuUtils.getPopupMenu(Arrays.asList(this), "Menu", te);
        menu.addSeparator();
        TrackMenuUtils.addDataRendererItems(menu, Arrays.asList("Heatmap", "Bar Chart", "Points", "Line Plot"), Arrays.asList(this));
        return menu;
    }

    /**
     * Return a value string for the tooltip window at the given location, or null to signal there is no value
     * at that location
     *
     * @param chr
     * @param position
     * @param mouseX
     * @param mouseY Mouse position relative to the enclosing panel
     * @param frame    @return
     */
    @Override
    public String getValueStringAt(String chr, double position, int mouseX, int mouseY, ReferenceFrame frame) {

        int trackY = mouseY - this.getY();
        int sampleIdx = trackY / sampleHeight;
        if (sampleIdx < 0 || sampleIdx >= sampleNames.size()) {
            return null;
        }
        String sampleName = sampleNames.get(sampleIdx);

        StringBuffer buf = new StringBuffer();
        LocusScore score = dataset.getSegmentAt(sampleName, chr, position, frame);
        // If there is no value here, return null to signal no popup
        if (score == null) {
            return null;
        }
        buf.append(sampleName + "<br>");
        if ((getDataRange() != null) && (getRenderer() instanceof XYPlotRenderer)) {
            buf.append("Data scale: " + getDataRange().getMinimum() + " - " + getDataRange().getMaximum() + "<br>");
        }

        buf.append(score.getValueString(position, mouseX, getWindowFunction()));
        return buf.toString();
    }

    /**
     * Get the score over the provided region for the given type. Different types
     * are processed differently. Results are cached according to the provided frameName,
     * if provided. If not, a string is created based on the inputs.
     *
     * @param chr
     * @param start
     * @param end
     * @param type
     * @return
     */
    @Override
    public void sortSamplesByValue(String chr, int start, int end, RegionScoreType type) {

        if (end <= start) {
            return;
        }

        // Compute a value for each sample
        Map<String, Float> sampleScores = new HashMap<>();
        int n = 0;
        for (String s : sampleNames) {
            List<LocusScore> scores = dataset.getSegments(s, chr);
            float regionScore = 0;
            int intervalSum = 0;
            boolean hasNan = false;
            for (LocusScore score : scores) {
                if ((score.getEnd() >= start) && (score.getStart() <= end)) {
                    int interval = Math.min(end, score.getEnd()) - Math.max(start, score.getStart());
                    float value = score.getScore();
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

        sampleNames.sort(new Comparator<String>() {
            @Override
            public int compare(String o1, String o2) {
                Float v1 = sampleScores.get(o1);
                Float v2 = sampleScores.get(o2);
                // Put NaNs at the end
                if (v1.isNaN()) {
                    return (v2.isNaN()) ? 0 : 1;
                } else if (v2.isNaN()) {
                    return -1;
                }
                return -Float.compare(v1, v2);
            }
        });
    }

    /**
     * The track itself is not filterable, but samples are.
     * @return
     */
    @Override
    public boolean isFilterable() {
        return false;
    }

    @Override
    public void filterSamples(TrackFilter trackFilter) {
        if(trackFilter == null || trackFilter.isShowAll())  {
            this.sampleNames = new ArrayList<>(dataset.getSampleNames());
        } else {
            this.sampleNames = trackFilter.evaluateSamples(dataset.getSampleNames());
        }
    }

}
