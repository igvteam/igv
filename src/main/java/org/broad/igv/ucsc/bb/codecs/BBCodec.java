package org.broad.igv.ucsc.bb.codecs;

import org.broad.igv.feature.BasicFeature;
import org.broad.igv.ucsc.bb.BedData;

public interface BBCodec {

    BasicFeature decode(BedData feat);
}
