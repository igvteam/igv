package org.broad.igv.bbfile.codecs;

import org.broad.igv.bbfile.BedData;
import org.broad.igv.feature.BasicFeature;

public interface BBCodec {

    BasicFeature decode(BedData feat);
}
