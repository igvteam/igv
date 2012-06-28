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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.renderer;

import org.apache.log4j.Logger;
import org.broad.igv.feature.IGVFeature;
import org.broad.tribble.Feature;

/**
 * @author jrobinso
 */
public abstract class FeatureRenderer implements Renderer<IGVFeature> {

    /**
     * Return the pixel position corresponding to the chromosomal position.
     */
    private static Logger log = Logger.getLogger(FeatureRenderer.class);

    private Feature highlightFeature = null;

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
