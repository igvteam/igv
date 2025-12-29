/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.renderer;

import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import htsjdk.tribble.Feature;

import java.awt.*;

/**
 * @author jrobinso
 */
public interface Renderer<T extends Feature> {

    void render(java.util.List<T> features, RenderContext context, Rectangle rect, Track track);

}
