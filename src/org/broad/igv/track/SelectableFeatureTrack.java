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

import org.broad.igv.feature.Exon;
import org.broad.igv.feature.IExon;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.renderer.SelectableFeatureRenderer;
import htsjdk.tribble.Feature;

import java.awt.event.MouseEvent;
import java.util.HashSet;
import java.util.Set;

/**
 * A FeatureTrack designed for selecting features/exons. Created to support Sashimi plot
 * window
 * User: jacob
 * Date: 2013-Jan-28
 */
public class SelectableFeatureTrack extends FeatureTrack {

    protected Set<IExon> selectedExons = new HashSet<IExon>();

    public SelectableFeatureTrack(FeatureTrack geneTrack) {
        super(geneTrack);
        this.setRendererClass(SelectableFeatureRenderer.class);
    }

    @Override
    public boolean handleDataClick(TrackClickEvent te) {
        MouseEvent e = te.getMouseEvent();
        Feature f = getFeatureAtMousePosition(te);

        //We allow any of these modifier keys for multi-select
        if (!e.isShiftDown() && !e.isControlDown() && !e.isMetaDown()) {
            clearSelectedExons();
        }

        boolean foundExon = false;
        if (f != null && f instanceof IGVFeature) {
            selectedFeature = (IGVFeature) f;
            double location = te.getFrame().getChromosomePosition(e.getX());
            if (selectedFeature.getExons() != null) {
                for (Exon exon : selectedFeature.getExons()) {
                    if (location >= exon.getStart() && location < exon.getEnd()) {
                        selectedExons.add(exon);
                        foundExon = true;
                        break;
                    }
                }
            }
        }

        ((SelectableFeatureRenderer) getRenderer()).setSelectedExons(selectedExons);

        return foundExon;
    }


    private void clearSelectedExons() {
        this.selectedExons.clear();
    }

    public Set<IExon> getSelectedExons() {
        return selectedExons;
    }
}

