package org.broad.igv.bedpe;

import java.io.IOException;
import java.util.List;

public interface InteractionSource {

    List<BedPE> getFeatures(String chr, int start, int end, double bpPerPixel, String normalization) throws IOException;

}