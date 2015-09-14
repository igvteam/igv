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

package org.broad.igv.renderer;

import com.google.common.base.Predicate;
import org.broad.igv.feature.Exon;
import org.broad.igv.feature.IExon;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.util.collections.CollUtils;

import java.awt.*;
import java.util.Collections;
import java.util.List;
import java.util.Set;

/**
 * User: jacob
 * Date: 2013-Jan-28
 */
public class SelectableFeatureRenderer extends IGVFeatureRenderer {

    private static final Stroke lineStroke;
    private static final Color selectedBorderColor;
    private Set<IExon> selectedExons = Collections.emptySet();



    static{
        lineStroke = new BasicStroke(2.0f);
        selectedBorderColor = Color.blue;
    }

    public SelectableFeatureRenderer(){
        AA_COLOR_1 = new Color(AA_COLOR_1.getRed(), AA_COLOR_1.getGreen(), AA_COLOR_1.getBlue(), 120);
        AA_COLOR_2 = new Color(AA_COLOR_2.getRed(), AA_COLOR_2.getGreen(), AA_COLOR_2.getBlue(), 120);
    }

    @Override
    public void render(List<IGVFeature> featureList, RenderContext context, Rectangle trackRectangle, Track track) {
        List<IGVFeature> featuresWithExons = CollUtils.filter(featureList, new Predicate<IGVFeature>() {
            @Override
            public boolean apply(IGVFeature input) {
                return input.getExons() != null && input.getExons().size() > 0;
            }
        });
        super.render(featuresWithExons, context, trackRectangle, track);
    }

    @Override
    protected void drawExonRect(Graphics blockGraphics, Exon exon, int x, int y, int width, int height) {
        boolean isSelected = selectedExons.contains(exon);
        for(IExon selEx: selectedExons){
            if(selEx.contains(exon.getStart()) || selEx.contains(exon.getEnd())){
                isSelected = true;
                break;
            }
        }

        if(isSelected){
            blockGraphics.clearRect(x, y, width, height);

            Graphics2D borderGraphics = (Graphics2D) blockGraphics.create();
            borderGraphics.setColor(selectedBorderColor);
            borderGraphics.setStroke(lineStroke);
            borderGraphics.drawRect(x, y, width, height);
        }else{
            super.drawExonRect(blockGraphics, exon, x, y, width, height);
        }
    }

    public void setSelectedExons(Set<IExon> selectedExons) {
        this.selectedExons = selectedExons;
    }
}
