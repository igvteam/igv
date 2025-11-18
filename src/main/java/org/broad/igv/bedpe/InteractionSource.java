package org.broad.igv.bedpe;

import java.io.IOException;
import java.util.List;

public interface InteractionSource {

    List<BedPE> getFeatures(String chr, int start, int end, double bpPerPixel, String normalization, int maxFeatureCount) throws IOException;

    default boolean hasNormalizationVector(String type, String chr, double bpPerPixel) {
        return false;
    }

    default List<String> getNormalizationTypes() {
        return null;
    }
}