/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.igv.renderer;

import org.igv.track.RenderContext;
import org.igv.track.Track;
import htsjdk.tribble.Feature;

import java.awt.*;

/**
 * @author jrobinso
 */
public interface Renderer<T extends Feature> {

    void render(java.util.List<T> features, RenderContext context, Rectangle rect, Track track);

}
