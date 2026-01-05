package org.igv.renderer;


import org.igv.logging.*;


/**
 * @author jrobinso
 */
public class GeneTrackRenderer extends IGVFeatureRenderer {
    Logger log = LogManager.getLogger(GeneTrackRenderer.class);

    /**
     * Method description
     *
     * @return
     */
    public String getDisplayName() {
        return "Genes";
    }
}
