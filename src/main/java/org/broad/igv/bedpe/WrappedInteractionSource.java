package org.broad.igv.bedpe;

import org.broad.igv.track.FeatureSource;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;

public class WrappedInteractionSource implements  InteractionSource {

    FeatureSource<BedPE> featureSource;


    public WrappedInteractionSource(FeatureSource<BedPE> featureSource) {
        this.featureSource = featureSource;
    }

    @Override
    public List<BedPE> getFeatures(String chr, int start, int end, double bpPerPixel, String normalization) throws IOException {
        Iterator<BedPE>  featureIterator = featureSource.getFeatures(chr, start, end);
        List<BedPE> list = new ArrayList<>();
        while (featureIterator != null && featureIterator.hasNext()) {
            list.add(featureIterator.next());
        }
        return list;
    }
}
