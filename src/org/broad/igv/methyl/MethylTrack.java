/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not
 * responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL), Version 2.1 which is
 * available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.methyl;

import org.broad.igv.Globals;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.renderer.PointsRenderer;
import org.broad.igv.track.AbstractTrack;
import org.broad.igv.track.RenderContext;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * Experimental class for specific dataset.
 * <p/>
 * if(locator.getPath().contains("RRBS_cpgMethylation") || locator.getPath().contains("")) {
 *
 * @author Jim Robinson
 * @date 4/19/12
 */
public class MethylTrack extends AbstractTrack {

    public static final int FIVE_MB = 5000000;
    public static final int FIFTY_MB = 50000000;

    private MethylDataSource dataSource;
    private Interval loadedInterval;
    private Renderer renderer;
    private int resolutionThreshold;
    private boolean loading = false;

    public MethylTrack(ResourceLocator dataResourceLocator, BBFileReader reader, Genome genome) throws IOException {
        super(dataResourceLocator);
        setHeight(60);
        renderer = new PointsRenderer();

        boolean isWGBS;
        if (reader.getAutoSql() != null && reader.getAutoSql().startsWith("table BisulfiteSeq")) {
            resolutionThreshold = FIVE_MB;
            dataSource = new CachingMethylSource(new BBMethylDataSource(reader, BBMethylDataSource.Type.USC, genome), resolutionThreshold);
            //dataSource = new BBMethylDataSource(reader, BBMethylDataSource.Type.USC, genome);
        } else {
            isWGBS = dataResourceLocator.getPath().contains("BiSeq_cpgMethylation");
            resolutionThreshold = isWGBS ? FIVE_MB : FIFTY_MB;
            dataSource = new CachingMethylSource(new BBMethylDataSource(reader, BBMethylDataSource.Type.ZILLER, genome), resolutionThreshold);
            //dataSource = new BBMethylDataSource(reader, BBMethylDataSource.Type.USC, genome);
        }

        loadedInterval = new Interval("", -1, -1, Collections.<MethylScore>emptyList());
        setDataRange(new DataRange(0, 100));
    }

    /**
     * Render the track in the supplied rectangle.  It is the responsibility of the track to draw within the
     * bounds of the rectangle.
     *
     * @param context the render context
     * @param rect    the track bounds, relative to the enclosing DataPanel bounds.
     */
    public void render(final RenderContext context, final Rectangle rect) {

        if (context.getChr().equals(Globals.CHR_ALL) || context.getScale() * rect.width > resolutionThreshold) {
            Graphics2D g = context.getGraphic2DForColor(Color.gray);
            Rectangle textRect = new Rectangle(rect);

            // Keep text near the top of the track rectangle
            textRect.height = Math.min(rect.height, 20);
            textRect.height = Math.min(rect.height, 20);
            String message = context.getChr().equals(Globals.CHR_ALL) ? "Zoom in to see features." :
                    "Zoom in to see features, or right-click to increase Feature Visibility Window.";
            GraphicUtils.drawCenteredText(message, textRect, g);
            return;
        }

        final String chr = context.getChr();
        final int start = (int) context.getOrigin();
        final int end = (int) context.getEndLocation();
        if (loadedInterval.contains(chr, start, end)) {
            renderer.render(loadedInterval.scores, context, rect, this);
        } else {
            if (!loading) {
                loading = true;
                Runnable runnable = new Runnable() {
                    public void run() {

                        try {
                            int width = (end - start) / 2;
                            int expandedStart = Math.max(0, start - width);
                            int expandedEnd = end + width;

                            List<MethylScore> scores = new ArrayList<MethylScore>(1000);
                            Iterator<MethylScore> iter = dataSource.query(chr, expandedStart, expandedEnd);
                            while (iter.hasNext()) {
                                scores.add(iter.next());
                            }
                            loadedInterval = new Interval(chr, expandedStart, expandedEnd, scores);
                            context.getPanel().repaint(); //rect);
                        } finally {
                            loading = false;
                        }
                    }
                };
                LongRunningTask.submit(runnable);
            }
        }
    }

    public Renderer getRenderer() {
        return renderer;  //To change body of implemented methods use File | Settings | File Templates.
    }


    static class Interval {
        String chr;
        int start;
        int end;
        List<MethylScore> scores;

        Interval(String chr, int start, int end, List<MethylScore> scores) {
            this.chr = chr;
            this.end = end;
            this.scores = scores;
            this.start = start;
        }

        boolean contains(String chr, int start, int end) {
            return this.chr.equals(chr) && this.start <= start && this.end >= end;
        }

    }
}
