package org.igv.bedpe;

import org.igv.track.RenderContext;

import java.util.List;

public interface BedPERenderer {
     void render(List<BedPE> features, RenderContext context, InteractionTrack.ArcOption arcOption);

}
