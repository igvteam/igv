package org.broad.igv.renderer;

/**
 * This class exists as a side effect of the fact that
 * we set renderers by class, rather than instantiating them.
 * If the latter were the case, we could just set the isGenotypeRenderer
 * parameter directly in a constructor or method.
 */
public class GenotypeRenderer extends IGVFeatureRenderer {

    public GenotypeRenderer() {
        super();
        isGenotypeRenderer = true;
        drawBoundary = false;
    }

}
