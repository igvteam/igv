package org.igv.track;

import org.igv.feature.Exon;
import org.igv.feature.IExon;
import org.igv.feature.IGVFeature;
import org.igv.renderer.SelectableFeatureRenderer;
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

    public SelectableFeatureTrack() {
        this.renderer = new SelectableFeatureRenderer();
    }

    public SelectableFeatureTrack(FeatureTrack geneTrack) {
        super(geneTrack);
        this.renderer = new SelectableFeatureRenderer();
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
            double location = te.getFrame().getChromosomePosition(e);
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

