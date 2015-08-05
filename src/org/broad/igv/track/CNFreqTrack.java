/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
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

package org.broad.igv.track;

import org.broad.igv.data.seg.FreqData;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.renderer.BarChartRenderer;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.util.Arrays;
import java.util.List;

/**
 * @author jrobinso
 * @date Oct 13, 2010
 */
public class CNFreqTrack extends AbstractTrack {


    FreqData data;
    BarChartRenderer renderer;

    public CNFreqTrack(ResourceLocator rl, String id, String name, FreqData fd) {
        super(rl, id, name);
        data = fd;

        float nSamples = data.getNumberOfSamples();
        this.setDataRange(new DataRange(-nSamples, 0, nSamples));
        this.setColor(Color.red);
        this.setAltColor(Color.blue);

        renderer = new BarChartRenderer();

        this.setMinimumHeight(25);
        this.setHeight(50);
        this.setSortable(false);

    }


    public void render(RenderContext context, Rectangle rect) {

        renderer.render(data.getDelCounts(context.getChr()), context, rect, this);
        renderer.render(data.getAmpCounts(context.getChr()), context, rect, this);
        renderer.setMarginFraction(0);
        renderer.renderBorder(this, context, rect);
        context.getGraphic2DForColor(Color.black).drawRect(rect.x, rect.y, rect.width, rect.height - 1);

    }


    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {

        List<LocusScore> ampScores = data.getAmpCounts(chr);
        List<LocusScore> delScores = data.getDelCounts(chr);
        StringBuffer buf = new StringBuffer();
        int startIdx = Math.max(0, FeatureUtils.getIndexBefore(position, ampScores));
        for (int i = startIdx; i < ampScores.size(); i++) {
            LocusScore ampScore = ampScores.get(i);
            if (position >= ampScore.getStart() && position <= ampScore.getEnd()) {
                buf.append("# of samples with log2(cn/2) &gt; &nbsp; " + data.getAmpThreshold() + ": ");
                buf.append(ampScore.getValueString(position, null));
                buf.append("<br># of samples with log2(cn/2) &lt;  " + data.getDelThreshold() + ":  ");
                buf.append(delScores.get(i).getValueString(position, null));
            }
        }
        return buf.length() == 0 ? null : buf.toString();
    }

    public Renderer getRenderer() {
        return renderer;
    }

    public boolean isLogNormalized() {
        return false;
    }


    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, String frameName) {
        return Integer.MIN_VALUE;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public IGVPopupMenu getPopupMenu(TrackClickEvent te) {

        IGVPopupMenu menu = new IGVPopupMenu();

        List<Track> selfAsList = Arrays.asList((Track) this);
        TrackMenuUtils.addSharedItems(menu, selfAsList, false);

        menu.addSeparator();
        menu.add(TrackMenuUtils.getRemoveMenuItem(selfAsList));


        return menu;
    }
}
