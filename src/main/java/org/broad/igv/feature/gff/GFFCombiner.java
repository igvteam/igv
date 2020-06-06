package org.broad.igv.feature.gff;

import htsjdk.tribble.Feature;
import org.broad.igv.feature.BasicFeature;

import java.util.Iterator;
import java.util.List;

public interface GFFCombiner {

    GFFCombiner addFeatures(Iterator<Feature> rawIter);

    void addFeature(BasicFeature bf);

    List<Feature> combineFeatures();
}
