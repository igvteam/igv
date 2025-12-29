/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.renderer;

import org.broad.igv.Globals;
import org.broad.igv.logging.*;
import org.broad.igv.feature.IGVFeature;
import htsjdk.tribble.Feature;

import java.awt.*;

/**
 * @author jrobinso
 */
public abstract class FeatureRenderer implements Renderer<IGVFeature> {

    /**
     * Return the pixel position corresponding to the chromosomal position.
     */
    private static Logger log = LogManager.getLogger(FeatureRenderer.class);

    protected boolean darkMode;
    protected final Color fontColor;

    private Feature highlightFeature = null;

    public FeatureRenderer() {
        this.darkMode = Globals.isDarkMode();
        this.fontColor = darkMode ? Color.WHITE : Color.BLACK;
    }

    public Feature getHighlightFeature() {
        return highlightFeature;
    }

    public void setHighlightFeature(Feature highlightFeature) {
        this.highlightFeature = highlightFeature;
    }

    /**
     * We render features 1 row at a time.
     * This is called once before rendering the N rows
     * Default implementation does nothing
     */
    public void reset() {
        //pass
    }
}
