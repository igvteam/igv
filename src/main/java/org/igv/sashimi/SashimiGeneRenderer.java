package org.igv.sashimi;

import org.igv.renderer.SelectableFeatureRenderer;

/**
 * Gene renderer for the Sashimi plot, positioning features through the plot's coordinate map.
 */
public class SashimiGeneRenderer extends SelectableFeatureRenderer {

    private final SashimiView view;

    public SashimiGeneRenderer(SashimiView view) {
        this.view = view;
    }

    @Override
    protected double getVirtualPixel(double chromosomeLocation, double origin, double locationScale) {
        return view.toPixel(chromosomeLocation);
    }
}
