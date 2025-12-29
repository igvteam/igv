/*
* To change this template, choose Tools | Templates
* and open the template in the editor.
*/
package org.igv.renderer;

//~--- non-JDK imports --------------------------------------------------------

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
