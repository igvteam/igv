package org.igv.seg;

import org.igv.feature.LocusScore;
import org.igv.track.DataType;
import org.igv.ui.panel.ReferenceFrame;

import java.util.List;

public interface SegmentedDataSource {

    DataType getType();

    double getDataMin();

    double getDataMax();

    List<String> getSampleNames();

    List<LocusScore> getFeatures(String sample, String chr);

    LocusScore getFeatureAt(String sample, String chr, double position, ReferenceFrame frame);
}

