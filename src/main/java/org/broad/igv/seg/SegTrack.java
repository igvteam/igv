package org.broad.igv.seg;

import htsjdk.tribble.Feature;
import org.broad.igv.data.DataSource;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.HeatmapRenderer;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.track.*;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.MouseableRegion;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class SegTrack extends AbstractTrack {

    SegmentedAsciiDataSet dataset;
    String sampleName;
    Color color = Color.RED;
    Color altColor = Color.BLUE;
    int sampleHeight = 15;
    List<String> sampleNames;
    Renderer renderer = new HeatmapRenderer();


    public SegTrack(ResourceLocator locator, String id, String name, SegmentedAsciiDataSet dataset, Genome genome) {
        super(locator, id, name);
        this.dataset = dataset;
        this.sampleNames = new ArrayList(dataset.getSampleNames());
        setDataRange(new DataRange(0, 1, 2));
        setTrackType(dataset.getType());
        initScale(dataset);
    }

    void initScale(SegmentedAsciiDataSet dataset) {

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

        int y = 0;
        for(String s : sampleNames) {
            this.sampleName = s;
            Rectangle r = new Rectangle(rect.x, rect.y + y, rect.width, sampleHeight);
            renderSample(context, r);
            y += sampleHeight;
        }
    }

    private void renderSample(RenderContext context, Rectangle rect) {
        List<LocusScore> scores = dataset.getSegments(sampleName, context.getChr());
        this.renderer.render(scores, context, rect, this);
    }

    @Override
    public void renderAttributes(Graphics2D graphics, Rectangle trackRectangle, Rectangle visibleRect, List<String> names, List<MouseableRegion> mouseRegions) {
        super.renderAttributes(graphics, trackRectangle, visibleRect, names, mouseRegions);
    }

    @Override
    public void renderName(Graphics2D g2D, Rectangle trackRectangle, Rectangle visibleRectangle) {
        super.renderName(g2D, trackRectangle, visibleRectangle);
    }

    @Override
    public boolean hasSamples() {
        return super.hasSamples();
    }

}
