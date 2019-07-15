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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.renderer;

import org.apache.log4j.Logger;
import org.broad.igv.feature.IGVFeature;
import htsjdk.tribble.Feature;

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
