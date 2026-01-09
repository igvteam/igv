package org.igv.bedpe;

import org.igv.track.FeatureSource;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;

public class WrappedInteractionSource implements InteractionSource {

    FeatureSource<BedPE> featureSource;


    public WrappedInteractionSource(FeatureSource<BedPE> featureSource) {
        this.featureSource = featureSource;
    }

    @Override
    public List<BedPE> getFeatures(String chr, int start, int end, double bpPerPixel, String normalization, int maxFeatureCount) throws IOException {
        Iterator<BedPE> featureIterator = featureSource.getFeatures(chr, start, end);
        List<BedPE> reservoir = new ArrayList<>();

        if (featureIterator == null) {
            return reservoir;
        }

        // If maxFeatureCount <= 0 treat as "no limit" and return all features
        if (maxFeatureCount <= 0) {
            while (featureIterator.hasNext()) {
                reservoir.add(featureIterator.next());
            }
            return reservoir;
        }

        int i = 0; // number of items seen so far

        while (featureIterator.hasNext()) {
            BedPE item = featureIterator.next();
            if (i < maxFeatureCount) {
                // fill the reservoir
                reservoir.add(item);
            } else {
                // replace elements with decreasing probability
                int r = ThreadLocalRandom.current().nextInt(i + 1); // 0..i
                if (r < maxFeatureCount) {
                    reservoir.set(r, item);
                }
            }
            i++;
        }

        return reservoir;
    }
}
