package org.broad.igv.bbfile.codecs;

import org.broad.igv.bbfile.BedFeature;
import org.broad.igv.feature.BasicFeature;

public interface BBCodec {

    BasicFeature decode(BedFeature feat);
}
