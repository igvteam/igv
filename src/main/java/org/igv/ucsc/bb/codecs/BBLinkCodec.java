package org.igv.ucsc.bb.codecs;

import org.igv.Globals;
import org.igv.feature.BasicFeature;
import org.igv.feature.IGVFeature;
import org.igv.ucsc.bb.BedData;

/**
 * Codec for UCSC "bigLink" files (bed4+1), the companion of bigChain.  Each record is one ungapped block of a chain,
 * named by the chain id.  Only the name is needed to join links to chains (see BBChainFeatureSource), so the record
 * is decoded to a minimal feature without the generic bed parsing.
 */
public class BBLinkCodec implements BBCodec {

    @Override
    public IGVFeature decode(BedData data) {
        BasicFeature feature = new BasicFeature(data.chr, data.chromStart, data.chromEnd);
        feature.setName(Globals.tabPattern.split(data.restOfFields, 2)[0]);
        return feature;
    }
}
