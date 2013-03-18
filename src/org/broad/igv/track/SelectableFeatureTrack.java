/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.track;

import org.broad.igv.feature.Exon;
import org.broad.igv.feature.IExon;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.renderer.SelectableFeatureRenderer;
import org.broad.tribble.Feature;

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

