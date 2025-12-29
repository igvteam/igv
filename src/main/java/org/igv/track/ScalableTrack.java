package org.igv.track;

import org.igv.feature.LocusScore;
import org.igv.renderer.DataRange;
import org.igv.ui.panel.ReferenceFrame;

import java.util.List;

/**
 * Created by jrobinson on 5/3/16.
 */
public interface ScalableTrack {

    Range getInViewRange(ReferenceFrame referenceFrame);

    void setDataRange(DataRange axisDefinition);
}
