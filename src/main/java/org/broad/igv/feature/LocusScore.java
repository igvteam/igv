/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.feature;

import org.broad.igv.track.WindowFunction;

/**
 * @author jrobinso
 */
public interface LocusScore extends htsjdk.tribble.Feature {

     float getScore();

    /**
     * Return a string to be used for popup text.   The WindowFunction is passed
     * in so it can be used to annotate the value.  The LocusScore object itself
     * does not "know" from what window function it was derived
     *
     * @param position       Zero-based genome position
     * @param mouseX
     * @param windowFunction
     * @return
     */
    default String getValueString(double position, int mouseX, WindowFunction windowFunction) {
        return "";
    }

}
