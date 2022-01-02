package org.broad.igv.ext.render;

import org.broad.igv.ext.IExtension;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.AlignmentBlock;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.sam.Gap;
import org.broad.igv.track.RenderContext;

public interface IIndelRenderingExtension extends IExtension {

    public void renderSmallInsertion(Alignment alignment,
                                     AlignmentBlock aBlock,
                                     RenderContext context,
                                     int h, int x, int y,
                                     AlignmentTrack.RenderOptions renderOptions);

    public void renderSmallInsertionWings(Alignment alignment,
                                          AlignmentBlock insertionBlock,
                                          RenderContext context,
                                          int pxH, int pxTop, int pxRight, int pxLeft,
                                          AlignmentTrack.RenderOptions renderOptions);

    public void renderDeletionGap(Alignment alignment,
                                  Gap gap,
                                  int y, int h, int x, int w,
                                  RenderContext context,
                                  AlignmentTrack.RenderOptions renderOptions);

}
