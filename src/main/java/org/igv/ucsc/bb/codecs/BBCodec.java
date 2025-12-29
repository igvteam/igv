package org.igv.ucsc.bb.codecs;

import org.igv.feature.BasicFeature;
import org.igv.feature.IGVFeature;
import org.igv.ucsc.bb.BedData;

public interface BBCodec {

    IGVFeature decode(BedData feat);
}
